#!/usr/bin/env dials.python

# 170209: Confirmed the result with prepare-cheetah-sacla-api2.

import dbpy
import stpy

import sys
import h5py
import math
import numpy as np
import re

XSIZE = 512
YSIZE = 1024
BL = 3

def write_crystfel_geom(filename, det_infos, energy):
    with open(filename, "w") as out:
        out.write("; CrystFEL geometry file produced by prepare_cheetah_api2.py\n")
        out.write(";   Takanori Nakane (takanori.nakane@bs.s.u-tokyo.ac.jp)\n")
        out.write("; for tiled but NOT reassembled images (512x8192 pixels)\n\n")
        out.write("; NOTE:\n")
        out.write("; This file is for multi-event HDF5 files. To use in single-event\n")
        out.write("; files, you have to edit 'data' record and prepare .beam file\n\n")
        out.write("clen = 0.0515  ; 51.5 mm camera length. You SHOULD optimize this!\n")
        out.write("res = 20000  ; 50 micron  1m /50 micron\n")
        out.write(";badrow_direction = x\n")
        out.write(";adu_per_eV = /LCLS/adu_per_eV\n")
        out.write(";max_adu = 250000 ; should NOT be used. see best practice on Web\n")
        out.write("data = /%/data ; change this to /LCLS/data for single-event files\n")
        out.write("photon_energy = /%%/photon_energy_ev ; roughly %f. change this to /LCLS/photon_energy_eV for single-event files\n\n" % energy)
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

	for i in range(8):
            det_info = det_infos[i]
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

            out.write("q%d/adu_per_eV = %f\n" % (det_id, 1.0 / (0.1 * energy))) # Keitaro's 0.1 photon
            out.write("q%d/min_fs = %d\n" % (det_id, 0))
            out.write("q%d/min_ss = %d\n" % (det_id, i * YSIZE))
            out.write("q%d/max_fs = %d\n" % (det_id, XSIZE - 1))
            out.write("q%d/max_ss = %d\n" % (det_id, (i + 1) * YSIZE - 1))
            out.write("q%d/fs = %fx %+fy\n" % (det_id, math.cos(rotation), math.sin(rotation)))
            out.write("q%d/ss = %fx %+fy\n" % (det_id, -math.sin(rotation), math.cos(rotation)))
            out.write("q%d/corner_x = %f\n" % (det_id, -detx))
            out.write("q%d/corner_y = %f\n\n" % (det_id, dety))

def write_cheetah_geom(filename, det_infos):
    posx = np.zeros((YSIZE * 8, XSIZE), dtype=np.float32)
    posy = posx.copy()
    posz = posx.copy()

    for i in range(8):
        det_info = det_infos[i]
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

def add_image(acc, readers, buffers, gains, tag_id, energy):
    for i in range(8):
        readers[i].collect(buffers[i], tag_id)
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

def run(runid):
    print "prepare-cheetah-sacla-api2.py version 1y0209"
    print " by Takanori Nakane at University of Tokyo"
    print

    # Get Run info
    run_info = dbpy.read_runinfo(BL, runid)
    high_tag = dbpy.read_hightagnumber(BL, runid)
    start_tag = run_info['start_tagnumber']
    end_tag = run_info['end_tagnumber']

    tag_list = dbpy.read_taglist_byrun(BL, runid)
    tag = tag_list[0]
    print "Run %d: HighTag %d, Tags %d (inclusive) to %d (exclusive), thus %d images" % (runid, high_tag, start_tag, end_tag, len(tag_list))
    print

    # Find detectors
    det_IDs = dbpy.read_detidlist(3, runid)
    print "Detector IDs: " + " ".join(det_IDs)
    det_IDs = sorted([x for x in det_IDs if re.match("MPCCD-8-.*-[1-8]", x)])
    print "MPCCD Octal IDs to use: " + " ".join(det_IDs)
    print

    # Get shutter status and find dark images
    shutter = [str2float(s) for s in dbpy.read_syncdatalist("xfel_bl_3_shutter_1_open_valid/status", high_tag, tag_list)]
    dark_tags = [tag for tag, is_open in zip(tag_list, shutter) if is_open == 0]
    print "Number of dark images to average: %d" % len(dark_tags)
    print

    # Setup buffer readers
    readers = [stpy.StorageReader(det_id, BL, (runid,)) for det_id in det_IDs]
    buffers = [stpy.StorageBuffer(reader) for reader in readers]
    
    # Read first image to get detector info
    for reader, buf in zip(readers, buffers):
        reader.collect(buf, dark_tags[0])
    det_infos = [buf.read_det_info(0) for buf in buffers]

    # Collect pulse energies
    pulse_energies  = [1000.0 * str2float(s) for s in dbpy.read_syncdatalist("xfel_bl_3_tc_spec_1/energy", high_tag, tuple(dark_tags))]
    print "xfel_bl_3_tc_spec_1/energy for %d frames:" % len(pulse_energies)
    print pulse_energies
    print
    mean_energy = np.mean(pulse_energies)
    print "Mean photon energy: %f eV" % mean_energy
    print "Configured photon energy: %f eV" % (1000.0 * dbpy.read_config_photonenergy(BL, runid))
    print

    # Create geometry files
    write_crystfel_geom("%d.geom" % runid, det_infos, mean_energy)
    write_cheetah_geom("%d-geom.h5" % runid, det_infos)

    # Create dark average
    print
    print "Calculating a dark average:"
    num_added = 0
    sum_buffer = np.zeros((YSIZE * 8, XSIZE), dtype=np.float64)
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
    sum_buffer[sum_buffer < 0] = 0
    sum_buffer[sum_buffer > np.iinfo(np.uint16).max] = np.iinfo(np.uint16).max
    averaged = sum_buffer.astype(np.uint16)

    f = h5py.File("%d-dark.h5" % runid, "w")
    f.create_dataset("data/data", data=averaged, compression="gzip", shuffle=True)
    f.close()
    print "Dark average was written to %s" % ("%d-dark.h5" % runid)

run(int(sys.argv[1]))
