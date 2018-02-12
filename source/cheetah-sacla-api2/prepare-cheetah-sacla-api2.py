#!/usr/bin/env dials.python

# 170209: Confirmed the result with prepare-cheetah-sacla-api2.

import dbpy
import stpy

import sys
import h5py
import math
import numpy as np
import re

VERSION = 180212
XSIZE = 512
YSIZE = 1024
NPANELS = 8

def log_error(message):
    with open("status.txt", "w") as out:
        out.write("Status: Status=Error-%s\n" % message)

def write_crystfel_geom(filename, det_infos, energy, clen):
    with open(filename, "w") as out:
        out.write("; CrystFEL geometry file produced by prepare_cheetah_api2.py\n")
        out.write(";   Takanori Nakane (takanori.nakane@bs.s.u-tokyo.ac.jp)\n")
        out.write("; for tiled but NOT reassembled images (512x8192 pixels)\n\n")
        out.write("clen = %.4f               ; %.1f mm camera length. You SHOULD optimize this!\n" % (clen * 1E-3, clen))
        out.write("res = 20000                 ; = 1 m /50 micron\n")
        out.write(";badrow_direction = x\n")
        out.write(";max_adu = 250000           ; should NOT be used. see best practice on CrystFEL's Web site\n")
        out.write("data = /%/data\n")
        out.write(";mask = /metadata/pixelmask ; this does not work in CrystFEL 0.6.2 (reported bug)\n")
        out.write("mask_good = 0x00            ; instead, we can specify bad regions below if necessary\n")
        out.write("mask_bad = 0xFF\n")
        out.write("photon_energy = /%%/photon_energy_ev ; roughly %.1f eV\n\n" % energy)
        out.write("; Definitions for geoptimiser\n")
        out.write("rigid_group_q1 = q1\n")
        out.write("rigid_group_q2 = q2\n")
        out.write("rigid_group_q3 = q3\n")
        out.write("rigid_group_q4 = q4\n")
        out.write("rigid_group_q5 = q5\n")
        out.write("rigid_group_q6 = q6\n")
        out.write("rigid_group_q7 = q7\n")
        out.write("rigid_group_q8 = q8\n\n")
        out.write("rigid_group_collection_connected = q1,q2,q3,q4,q5,q6,q7,q8\n")
        out.write("rigid_group_collection_independent = q1,q2,q3,q4,q5,q6,q7,q8\n\n")

        out.write("; Panel definitions\n")
        for i, det_info in enumerate(det_infos):
            gain = det_info['mp_absgain']
            detx = det_info['mp_posx']
            dety = det_info['mp_posy']
            detz = det_info['mp_posz']
            rotation = det_info['mp_rotationangle'] * (math.pi / 180.0) # rad
            pixel_size = det_info['mp_pixelsizex']
            print "gain %f pos (%f, %f, %f) rotation %f energy %f" % (gain, detx, dety, detz, rotation, energy)

            detx /= pixel_size; dety /= pixel_size;
            det_id = i + 1

            # Nphotons = S [ADU] * G [e-/ADU] / (E [eV] * 3.65 [eV/e-]) according to the manual.
            # Thus, ADU/eV = 1/(3.65*G)
            out.write("; sensor %s\n" % det_info['id'])
            out.write("q%d/adu_per_eV = %f\n" % (det_id, 1.0 / (0.1 * energy))) # Keitaro's 0.1 photon
            out.write("q%d/min_fs = %d\n" % (det_id, 0))
            out.write("q%d/min_ss = %d\n" % (det_id, i * YSIZE))
            out.write("q%d/max_fs = %d\n" % (det_id, XSIZE - 1))
            out.write("q%d/max_ss = %d\n" % (det_id, (i + 1) * YSIZE - 1))
            out.write("q%d/fs = %fx %+fy\n" % (det_id, math.cos(rotation), math.sin(rotation)))
            out.write("q%d/ss = %fx %+fy\n" % (det_id, -math.sin(rotation), math.cos(rotation)))
            out.write("q%d/corner_x = %f\n" % (det_id, -detx))
            out.write("q%d/corner_y = %f\n\n" % (det_id, dety))

        border, outer_border = get_border(det_infos[0]['id'])
        if border != 0:
            out.write("; Bad regions near edges of each sensor\n")
            # NOTE: ranges are inclusive in CrystFEL
            # left
            out.write("badv1/min_fs = %d\n"   % 0)
            out.write("badv1/max_fs = %d\n"   % (border - 1))
            out.write("badv1/min_ss = %d\n"   % 0)
            out.write("badv1/max_ss = %d\n\n" % (YSIZE * NPANELS - 1))
            # right
            out.write("badv2/min_fs = %d\n"   % (XSIZE - border))
            out.write("badv2/max_fs = %d\n"   % (XSIZE - 1))
            out.write("badv2/min_ss = %d\n"   % 0)
            out.write("badv2/max_ss = %d\n\n" % (YSIZE * NPANELS - 1))
            for i in xrange(NPANELS):
                out.write("badq%dh1/min_fs = %d\n"   % (i, 0))
                out.write("badq%dh1/max_fs = %d\n"   % (i, XSIZE - 1))
                out.write("badq%dh1/min_ss = %d\n"   % (i, YSIZE * i))
                out.write("badq%dh1/max_ss = %d\n\n" % (i, YSIZE * i + border - 1))
        if outer_border != 0:
            out.write("; Bad regions near outer edges of each sensor due to amplifier shields\n")
            out.write(";  you might want to optimize these widths (edit min_ss)\n")
            for i in xrange(NPANELS):
                out.write("badq%dh2/min_fs = %d\n"   % (i, 0))
                out.write("badq%dh2/max_fs = %d\n"   % (i, XSIZE - 1))
                out.write("badq%dh2/min_ss = %d\n"   % (i, YSIZE * (i + 1) - outer_border))
                out.write("badq%dh2/max_ss = %d\n\n" % (i, YSIZE * (i + 1) - 1))

def write_cheetah_geom(filename, det_infos):
    posx = np.zeros((YSIZE * NPANELS, XSIZE), dtype=np.float32)
    posy = posx.copy()
    posz = posx.copy()

    for i, det_info in enumerate(det_infos):
        gain = det_info['mp_absgain']
        detx = det_info['mp_posx'] * 1E-6 # m
        dety = det_info['mp_posy'] * 1E-6
        detz = det_info['mp_posz'] * 1E-6
        rotation = det_info['mp_rotationangle'] * (math.pi / 180.0) # rad
        pixel_size = det_info['mp_pixelsizex'] * 1E-6 # m
        
        fast_x = math.cos(rotation) * pixel_size
        fast_y = math.sin(rotation) * pixel_size;
        slow_x = -math.sin(rotation) * pixel_size
        slow_y = math.cos(rotation) * pixel_size;

        y = np.arange(YSIZE)
        for x in xrange(XSIZE):
            posx[i * YSIZE + y, x] = fast_x * x + slow_x * y - detx
            posy[i * YSIZE + y, x] = -(fast_y * x + slow_y * y + dety)
            posz[i * YSIZE + y, x] = 0

    f = h5py.File(filename, "w")
    f.create_dataset("x", data=posx, compression="gzip", shuffle=True)
    f.create_dataset("y", data=posy, compression="gzip", shuffle=True)
    f.create_dataset("z", data=posz, compression="gzip", shuffle=True)
    f.close()

def get_border(det_name):
    if re.match("MPCCD-8B", det_name): # Phase 3 detector
        return (5, 23) # based on 17Jul-P3Lys @ 10 keV
    if re.match("MPCCD-8N", det_name): # Compact detector with amp shields
        return (0, 22) # based on 17Jul-Kuma @ 7 keV
    else:
        return (0, 0)

def make_pixelmask(borders):
    border, outer_border = borders

    mask = np.zeros((YSIZE * NPANELS, XSIZE), dtype=np.uint16)
    mask[:, 0:border] = 1
    mask[:, (XSIZE - border):XSIZE] = 1

    for i in xrange(NPANELS):
        mask[(YSIZE * i):(YSIZE * i + border), :] = 1
        mask[(YSIZE * (i + 1) - outer_border):(YSIZE * (i + 1)), :] = 1

    return mask

def write_metadata(filename, det_infos, clen, comment):
    f = h5py.File(filename, "w")
    
    f["/metadata/pipeline_version"] = VERSION
    f["/metadata/run_comment"] = comment
    f["/metadata/sensor_id"] = [det_info['id'] for det_info in det_infos]
    f["/metadata/posx_in_um"] = [det_info['mp_posx'] for det_info in det_infos]
    f["/metadata/posy_in_um"] = [det_info['mp_posy'] for det_info in det_infos]
    f["/metadata/posz_in_um"] = [det_info['mp_posz'] for det_info in det_infos]
    f["/metadata/angle_in_rad"] = [det_info['mp_rotationangle'] for det_info in det_infos]
    f["/metadata/pixelsizex_in_um"] = [det_info['mp_pixelsizex'] for det_info in det_infos]
    f["/metadata/pixelsizey_in_um"] = [det_info['mp_pixelsizey'] for det_info in det_infos]
    f["/metadata/distance_in_mm"] = clen
    pixel_mask = make_pixelmask(get_border(det_infos[0]['id']))
    f.create_dataset("/metadata/pixelmask", data=pixel_mask, compression="gzip", shuffle=True)
    f.close()

def add_image(acc, readers, buffers, gains, tag_id, energy):
    for i in range(NPANELS):
        try:
            readers[i].collect(buffers[i], tag_id)
        except:
            log_error("FailedOn_collect_data")
            sys.exit(-1)
        data = buffers[i].read_det_data(0)
        gain = gains[i] * 3.65 / 0.1 / energy
        acc[(YSIZE * i):(YSIZE * (i + 1)),] += data * gain

    return 1

def str2float(str):
    m = re.match("-?\d+(.\d+)?(e[+-]?\d+)?", str)
    if m is not None:
        return float(m.group(0))
    else:
        return None

def run(runid, bl=3, clen=50.0):
    # Beamline specific constants
    if bl == 2:
        sensor_spec = "xfel_bl_2_tc_spec_1/energy"
        sensor_shutter = "xfel_bl_2_shutter_1_open_valid/status"
    elif bl == 3:
        sensor_spec = "xfel_bl_3_tc_spec_1/energy"
        sensor_shutter = "xfel_bl_3_shutter_1_open_valid/status"
    else:
        log_error("BadBeamline")
        sys.exit(-1)

    # Get Run info
    try:
        run_info = dbpy.read_runinfo(bl, runid)
    except:
        log_error("BadRunID")
        sys.exit(-1)
    high_tag = dbpy.read_hightagnumber(bl, runid)
    start_tag = run_info['start_tagnumber']
    end_tag = run_info['end_tagnumber']

    tag_list = dbpy.read_taglist_byrun(bl, runid)
    tag = tag_list[0]
    print "Run %d: HighTag %d, Tags %d (inclusive) to %d (exclusive), thus %d images" % (runid, high_tag, start_tag, end_tag, len(tag_list))
    comment = dbpy.read_comment(bl, runid)
    print "Comment: %s" % comment
    print

    # Find detectors
    det_IDs = dbpy.read_detidlist(bl, runid)
    print "Detector IDs: " + " ".join(det_IDs)
    det_IDs = sorted([x for x in det_IDs if re.match("^MPCCD-8.*-[1-8]$", x)])
    if len(det_IDs) != 8:
        log_error("NoSupportedDetectorFound")
        sys.exit(-1)
    print "MPCCD Octal IDs to use: " + " ".join(det_IDs)
    print

    # Get shutter status and find dark images
    try:
        shutter = [str2float(s) for s in dbpy.read_syncdatalist(sensor_shutter, high_tag, tag_list)]
    except:
        log_error("NoShutterStatus")
        sys.exit(-1)
    dark_tags = [tag for tag, is_open in zip(tag_list, shutter) if is_open == 0]
    
    if bl == 2 and runid >= 32348 and runid <= 33416:
	# 2018 Feb: Unreliable shutter status. Should use PD and take darks only at the beginning of a run
        print "The shutter status was unreliable for runs in 2018 Feb."
        print "The number of tags with shutter closed:", len(dark_tags)
        print "Since the above value is not reliable, we use X-ray PD values instead."
        xray_pd = "xfel_bl_2_st_3_bm_1_pd/charge"
        pd_values = [str2float(s) for s in dbpy.read_syncdatalist(xray_pd, high_tag, tag_list)]
        dark_tags = []
        is_head = True
        for tag, pd in zip(tag_list, pd_values):
             if pd == None and is_head:
                 dark_tags.append(tag)
             else:
                 is_head = False
        print "Number of tags without X-ray:", len([1 for pd_val in pd_values if pd_val is None])
        print "But we use only tags at the beginning of a run."

    if len(dark_tags) == 0:
        log_error("NoDarkImage")
        sys.exit(-1)
    print "Number of dark images to average: %d" % len(dark_tags)
    print

    # Setup buffer readers
    try:
        readers = [stpy.StorageReader(det_id, bl, (runid,)) for det_id in det_IDs]
    except:
        log_error("FailedOn_create_streader")
        sys.exit(-1)
    try:
        buffers = [stpy.StorageBuffer(reader) for reader in readers]
    except:
        log_error("FailedOn_create_stbuf")
        sys.exit(-1)
    
    # Read first image to get detector info
    det_infos = []
    for reader, buf in zip(readers, buffers):
        try:
            reader.collect(buf, dark_tags[0])
        except:
            log_error("FailedOn_collect_data")
            sys.exit(-1)
    det_infos = [buf.read_det_info(0) for buf in buffers]
    for i, det_info in enumerate(det_infos):
        det_info['id'] = det_IDs[i]

    # Collect pulse energies
    config_photon_energy = 1000.0 * dbpy.read_config_photonenergy(bl, runid)
    config_photon_energy_sensible = True
    if config_photon_energy < 5000 or config_photon_energy > 14000:
        print "WARNING: dbpy.read_config_photonenergy returned %f eV, which is absurd!" % config_photon_energy
        print "         Report this to SACLA DAQ team."
        print "         This is not problematic unless the inline spectrometer is also broken." 
        config_photon_energy_sensible = False

    pulse_energies_in_keV  = [str2float(s) for s in dbpy.read_syncdatalist(sensor_spec, high_tag, tuple(dark_tags))]
    pulse_energies = []
    for tag, energy in zip(dark_tags, pulse_energies_in_keV):
        if energy is not None and energy > 0:
            pulse_energies.append(energy * 1000.0)
        else:
            print "WARNING: The wavelength from the inline spectrometer does not look sensible for tag %d." % tag
            if config_photon_energy_sensible:
                pulse_energies.append(config_photon_energy)
                print "         Used the accelerator config value instead."
            else:
                pulse_energies.append(7000.0)
                print "         The accelerator config value is also broken; assumed 7 keV as a last resort!"

    print
    mean_energy = np.mean(pulse_energies)
    print "Mean photon energy: %f eV" % mean_energy
    print "Configured photon energy: %f eV" % config_photon_energy
    print

    # Create geometry files
    write_crystfel_geom("%d.geom" % runid, det_infos, mean_energy, clen)
    write_cheetah_geom("%d-geom.h5" % runid, det_infos)

    # Write metadata
    write_metadata("%d.h5" % runid, det_infos, clen, comment)

    # Create dark average
    print
    print "Calculating a dark average:"
    num_added = 0
    sum_buffer = np.zeros((YSIZE * NPANELS, XSIZE), dtype=np.float64)
    gains = [det_info['mp_absgain'] for det_info in det_infos]

    for j, tag_id in enumerate(dark_tags):
        print "Processing tag %d (%2.1f%% done)" % (tag_id, 100.0 * (j + 1) / len(dark_tags))
        if (j % 5 == 0):
            with open("status.txt", "w") as status:
                status.write("Status: Total=%d,Processed=%d,Status=DarkAveraging\n" % (len(dark_tags), j + 1))
        num_added += add_image(sum_buffer, readers, buffers, gains, tag_id, pulse_energies[j])
    print "\nDone. Averaged %d frames." % num_added
  
    if (num_added < 1):
        return -1

    sum_buffer /= num_added
    ushort_max = np.iinfo(np.uint16).max
    print " #neg (< 0) %d, #overflow (> %d) %d" % (np.sum(sum_buffer < 0), ushort_max, np.sum(sum_buffer > ushort_max))
    sum_buffer[sum_buffer < 0] = 0
    sum_buffer[sum_buffer > ushort_max] = ushort_max
    averaged = sum_buffer.astype(np.uint16)

    # In the Phase 3 detector, some pixels average to negative values.
    # Most are around -0.1 and all those below -1 are at panel edges that will be masked.
    # So we don't have to worry about them.

    f = h5py.File("%d-dark.h5" % runid, "w")
    f.create_dataset("/data/data", data=averaged, compression="gzip", shuffle=True)
#    f.create_dataset("/data/raw", data=sum_buffer, compression="gzip", shuffle=True)
    f.close()
    print "Dark average was written to %s" % ("%d-dark.h5" % runid)

import optparse
parser = optparse.OptionParser()

parser.add_option("--bl", dest="bl", type=int, default=3, help="Beamline")
parser.add_option("--clen", dest="clen", type=float, default=50.0, help="Camera distance")
opts, args = parser.parse_args()

if (opts.bl != 2 and opts.bl !=3):
    print "--bl must be 2 or 3."
    sys.exit(-1)

if len(args) != 1:
    print "Usage: prepare-cheetah-sacla-api2.py runid [--bl 3] [--clen 50.0]"
    sys.exit(-1)
runid = int(args[0])

print "prepare-cheetah-sacla-api2.py version %d" % VERSION
print " by Takanori Nakane at University of Tokyo"
print
print "Option: bl               = %d" % opts.bl
print "Option: clen             = %.1f mm" % opts.clen
print

run(runid=runid, bl=opts.bl, clen=opts.clen)
