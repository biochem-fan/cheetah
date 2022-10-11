# Ref:
# http://stackoverflow.com/questions/21579658/embedding-a-live-updating-matplotlib-graph-in-wxpython
# http://stackoverflow.com/questions/4098131/how-to-update-a-plot-in-matplotlib

import re
import sys
import wx
import numpy as np
import matplotlib.figure as mfigure
import matplotlib.animation as manim
import optparse

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg

class RingBuffer():
    def __init__(self, nmax):
        self.nmax = nmax
        self.buf = np.zeros(self.nmax)
        self.n = 0

    def append(self, x):
        if self.n == self.nmax:
            self.buf = np.roll(self.buf, -1)
        else:
            self.n += 1

        self.buf[self.n - 1] = x

    def clear(self):
        self.n = 0
        self.buf.fill(0)

    def get(self, return_all=False):
        head = self.buf[0:self.n]

        if return_all:
            if self.n == 0:
                return np.zeros(self.nmax)
            else:
                return np.concatenate((head, np.zeros(self.nmax - self.n)))
        else:
            return head

    def __len__(self):
        return self.n

class PlotWindow(wx.Frame):
    TAGS_IN_A_RUN = 5151

    def __init__(self, parent, opts):
        super(PlotWindow,self).__init__(None, wx.ID_ANY, size=(1024, 550), title="Hitrate Plotter")

        self.hit_window = opts.window
        self.hit_threshold = opts.threshold
        self.fixed_filename = opts.filename
        self.filename = self.fixed_filename
        self.saturation = opts.saturation

        self.runid = 355800

        self.fig = mfigure.Figure()
        self.fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.95, wspace=0.15, hspace=0.15)
        self.ax = self.fig.add_subplot(111)
        self.ax_hit = self.ax.twinx()
        self.canvas = FigureCanvasWxAgg(self, wx.ID_ANY, self.fig)

        # Data
        self.peaks = RingBuffer(self.TAGS_IN_A_RUN)
        self.framenumber = RingBuffer(self.TAGS_IN_A_RUN)
        self.saturated = RingBuffer(self.TAGS_IN_A_RUN)
        self.frame_sat = RingBuffer(self.TAGS_IN_A_RUN)

        # Parser internals
        self.frame_read_bytes = 0
        self.peak_read_bytes = 0
        self.peak_cur_tag = -1
        self.peak_nsat = 0

        self.ax.set_ylim([0,300])
        self.ax_hit.set_ylim([0, 100])
        self.ax.set_xlabel("Tag ID")
        self.ax.set_ylabel("Number of peaks")
        self.ax_hit.set_ylabel("Hit rate (%)")
        self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.ax.get_xaxis().get_major_formatter().set_scientific(False) # no exponential notation

        self.peakdots, = self.ax.plot(0, 0, color="r", marker="o", linestyle="None", markersize=3, markeredgewidth=0)
        self.satdots, = self.ax.plot(0, 0, color="b", marker="o", linestyle="None", markersize=3, markeredgewidth=0)
        self.hitline, = self.ax_hit.plot(0, 0, 'b-')
        self.animator = manim.FuncAnimation(self.fig,self.anim, interval=3000) # 3 sec. cf. saveInterval in cheetah

        self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
#       self.text_runid = wx.TextCtrl(self)
#       self.label_runid = wx.StaticText(self, wx.ID_ANY, "XXXXXX")
        self.back_button = wx.Button(self, wx.ID_ANY, " < ")
        self.back_button.Bind(wx.EVT_BUTTON, self.onBackPush)
        self.next_button = wx.Button(self, wx.ID_ANY, " > ")
        self.next_button.Bind(wx.EVT_BUTTON, self.onNextPush)

        self.hsizer.Add(self.back_button, 0, wx.ALIGN_CENTER_VERTICAL)        
#        self.hsizer.Add(self.label_runid, 0, wx.ALIGN_CENTER_VERTICAL)
        self.hsizer.Add(self.next_button, 0, wx.ALIGN_CENTER_VERTICAL)
#        self.hsizer.Add(self.text_runid, 1, wx.ALL, 5)
#        self.next_button.Enable(False)
#        self.back_button.Enable(False)

        self.vsizer = wx.BoxSizer(wx.VERTICAL)
        self.vsizer.Add(self.canvas, 0, wx.EXPAND | wx.RIGHT | wx.TOP)
        self.vsizer.Add(self.hsizer, 0, wx.EXPAND | wx.RIGHT)
        self.SetSizer(self.vsizer)

        self.Show()

    def onBackPush(self, evt):
        self.runid -= 1
        self.anim(0)

    def onNextPush(self, evt):
        self.runid += 1
        self.anim(0)

    def getLatest(self):
        import os
        import glob

        frames = glob.glob("frames-run*.txt")
        frames.sort(key=os.path.getmtime, reverse=True)
        return frames[0]

    def anim(self, i):
#        if True:
#            filename = "frames-run%6d.txt" % self.runid
        if self.fixed_filename is None:
            new_filename = self.getLatest()
            if self.filename != new_filename:
                self.frame_read_bytes = 0
                self.peak_read_bytes = 0
                self.filename = new_filename
                
        m = re.search("run(\d+)", self.filename)
        if m is not None:
            self.runid = int(m.group(1))
       
        self.SetTitle(self.filename)

        self.parse_frames(self.filename)
        if len(self.peaks) == 0:
            return
        peaks = self.peaks.get()
        framenumber = self.framenumber.get()
        order = framenumber.argsort()
        peaks = peaks[order]
        framenumber = framenumber[order]

        self.ax.set_ylim([0, np.max(peaks)])
        self.ax.set_xlim([framenumber[0], framenumber[0] + self.TAGS_IN_A_RUN * 2])

        self.peakdots.set_xdata(framenumber)
        self.peakdots.set_ydata(peaks)

        hitrate = np.convolve(100 * (peaks[0:len(peaks)] > self.hit_threshold), 
                              np.ones(self.hit_window) / self.hit_window, mode="valid")
        self.hitline.set_xdata(framenumber[0:len(hitrate)])
        self.hitline.set_ydata(hitrate)

        self.parse_peaks(self.filename.replace("frames", "peaks"), self.saturation)
        self.satdots.set_xdata(self.frame_sat.get())
        self.satdots.set_ydata(self.saturated.get())

        self.fig.canvas.draw()
 
    def parse_frames(self, filename):
        try:
            import os
            size = os.path.getsize(filename)
        except:
            size = 0
        # This is necessary to limit latency on some file systems.

        cnt = 0
        with open(filename, "r") as f:
            f.seek(self.frame_read_bytes)

            for line in f:
                cnt += 1
                self.frame_read_bytes += len(line)
                
                if line[0] == "#":
                    continue
                columns = line.split(",")
                if len(columns) != 19:
                    continue

                try:
                    framenumber, peaks = int(columns[1]), int(columns[11])
                    self.framenumber.append(framenumber)
                    self.peaks.append(peaks)
                except:
                    print("Parse error in %s line: %s" % (filename, line))

                if self.frame_read_bytes > size:
                    break

    def parse_peaks(self, filename, saturation=2000):
        try:
            import os
            size = os.path.getsize(filename)
        except:
            size = 0

        cnt = 0
        with open(filename, "r") as f:
            f.seek(self.peak_read_bytes)

            for line in f:
                cnt += 1
                self.peak_read_bytes += len(line)

                if line[0] == "#":
                    continue
                columns = line.split(",")
                if len(columns) < 16:
                    continue

                try:
                        tag, maxi = int(columns[0]), float(columns[13])
        	        
                        if self.peak_cur_tag != tag:
                            if self.peak_cur_tag != -1:
                                self.frame_sat.append(self.peak_cur_tag)
                                self.saturated.append(self.peak_nsat)
                            self.peak_nsat = 0
                            self.peak_cur_tag = tag
                        if maxi >= saturation:
                            self.peak_nsat += 1
                except:
                    print("Parse error in %s line: %s" % (filename, line))

                if self.peak_read_bytes > size:
                    break

            # FIXME: This ignores the last frame?

print()
print("Cheetah Online plotter version 2022/10/11")
print("   by Takanori Nakane")
print()

parser = optparse.OptionParser()
parser.add_option("--window", dest="window", type=int, default=30, help="window size for moving average")
parser.add_option("--threshold", dest="threshold", type=int, default=20, help="minimum number of peaks for hits")
parser.add_option("--filename", dest="filename", type=str, default=None, help="filename to display (set None to follow current directory)")
parser.add_option("--saturation", dest="saturation", type=int, default=15000, help="saturation threshold")
opts, args = parser.parse_args()

print("Option: window = %d" % opts.threshold)
print("Option: threshold = %d" % opts.window)
print("Option: filename = %s" % opts.filename)
print("Option: saturation = %d" % opts.saturation)

app = wx.App(False)
frame = PlotWindow(None, opts)
app.MainLoop()
