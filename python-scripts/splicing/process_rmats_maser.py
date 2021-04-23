import argparse
import os
import glob
from collections import defaultdict
import re
import numpy as np
import math
import matplotlib.pyplot as plt
from typing import Union


def read_groups(groups: str) -> dict:
    """
    Gather differential splicing files from
    several rmats runs. It looks for files 
    provided in the directories listed in 
    the `groups` argument.

    :param str groups: File mapping a run 
    directory (1st col) to the corresponding 
    runID/comparison (e.g. ctrl_vs_treatment)

    :return dict: Dictionary grouping a comparison
    name with corresponding rmats output files
    """
    events = ["SE", "MXE", "A3SS", "A5SS", "RI"]
    d = {}
    with open(groups, 'r') as f:
        for line in f:
            _l = line.rstrip()
            items = _l.split("\t")
            if os.path.isdir(items[0]):
                _dir, group = items[0], items[1]
                d[_dir] = group
            else:
                raise ValueError("{} directory does not exist.".format(items[0]))

    groups_dict = {}
    for _dir, group in d.items():
        files_in_dir = []
        for ev in events:
            files_in_dir.extend(glob.glob(os.path.join(_dir, "*sign*{}*".format(ev))))

        groups_dict[group] = files_in_dir

    return groups_dict


def get_clean_coordinates(coord, event_type):
    """
    Processes coordinates of the event
    :param str coord: Input coordinate
    :param str event_type: Type of event

    :return str: Target coordinates for a single event.
    For a SE, it's the cassette exon, ALTSS events, it is
    the long exon form, for IR, it's the exon with the intron
    retained, and for the MXE, it's the coordinates of the
    two mutually exclusive exons.

    :return str: Spanning coordinates for a single event
    """
    coord_fields = coord.split(":")
    chrom = coord_fields[0]
  
    coord_fields = [item for sublist in [elem.split("-") for elem in coord_fields[2:]] for item in sublist]

    if event_type == "MXE":
        target_exon = "{}:{}-{};{}:{}-{}".format(chrom, coord_fields[0], coord_fields[1],
                                                 chrom, coord_fields[2], coord_fields[3])
    
    else:
        target_exon = "{}:{}-{}".format(chrom, coord_fields[0], coord_fields[1])
    spanning_coord = "{}:{}-{}".format(chrom, min(coord_fields), max(coord_fields))
    
    return target_exon, spanning_coord


def get_relevant_fields(file: str,
                        event_type: str,
                        sign_events: defaultdict, 
                        psi_threshold: float,
                        fdr_threshold: float,
                        group_name: str = None):
    """
    Processes individual event files from a
    single rMATS run.

    :param str file: Input file to be analyzed
    :param str event_type: Event type of the
    file to be analyzed.
    :param defaultdict sign_events: Dictionary
    to update significant events per comparison
    :param float psi_threshold: PSI thresh.
    :param float fdr_threshold: FDR thresh.
    :param str group_name: Name of the group
    comparison, if any.

    :return defaultdict: Updated dict
    """
    event_type_map = {"SE": "ES"}
   
    with open(file, 'r') as infile:
        next(infile)
        for line in infile:
            l = line.rstrip().split("\t")
            geneid = l[1].split(".")[0]
            genesymbol = l[2]
            pvalue = l[3]
            fdr = l[4]

            # rmats by default display the variation of
            # group1 vs group2, as opposed to the
            # other methods
            deltapsi = str(-float(l[5]))
            psi_1 = l[6]
            psi_2 = l[7]
            coordinates = ':'.join(l[8:14])
            if abs(float(deltapsi)) < psi_threshold or float(fdr) > fdr_threshold:
                continue

            target_exon, spanning_coord = get_clean_coordinates(coordinates, event_type)
            ev = event_type_map.get(event_type, event_type)
            if group_name:
                sign_events[coordinates].append([genesymbol, geneid, ev, deltapsi, pvalue, fdr, psi_1, psi_2,
                                                 target_exon, spanning_coord, group_name])
            else:
                sign_events[coordinates].append([genesymbol, geneid, ev, deltapsi, pvalue, fdr, psi_1, psi_2,
                                                 target_exon, spanning_coord])
    infile.close()
    return sign_events


def process_rmats_files(input_data: Union[str, dict],
                        psi_threshold: float,
                        fdr_threshold: float):
    """
    Process rMATS files so that readable tables are output

    :param Union[str, dict] input_data: Input data from
    rMATS runs. If a str is given (path where rMATS files
    are located), it means that a single is to be analyzed.
    If a dictionary is given instead, it means that multiple
    comparison are going to be analyzed, so the files for
    each comparison were already identified by reading groups before.

    :param float psi_threshold: PSI threshold to consider an
    event as differentially spliced
    :param float fdr_threshold: FDR threshold to consider an
    event as differentially spliced

    :return:
    """
    events = ["SE", "MXE", "A3SS", "A5SS", "RI"]
    header = ["#coordinates_event_id", "gene_name", "gene_id", "event_type", "deltaPSI", "pvalue", "FDR",
              "PSI_1", "PSI_2", "target_coordinates", "spanning_coordinates"]

    sign_events = defaultdict(list)

    if isinstance(input_data, dict):
        header.append("group")

        for group, files in input_data.items():
            for f in files:
                event_type = [ev for ev in events if ev in f]
                if len(event_type) > 1:
                    raise ValueError('Multiple event type flags found in the {} file ({}).'.format(f, event_type))
                elif len(event_type) == 0:
                    raise ValueError('No rMATS event type tag found in the {} file.'.format(f))
                else:
                    event_type = event_type[0]
                    sign_events = get_relevant_fields(f, event_type,
                                                      sign_events, psi_threshold,
                                                      fdr_threshold, group_name=group)
    else:
        for ev in events:
            file = glob.glob(os.path.join(input_data, "*sign*{}*".format(ev)))
            if len(file) == 1:
                f = file[0]
                sign_events = get_relevant_fields(f, ev, sign_events,
                                                  psi_threshold, fdr_threshold)
            elif len(file) > 1:
                raise ValueError('Multiple {} files found in the input directory {}'.format(ev, input_data))

    return sign_events, header


def write_output(events: defaultdict, header: str, 
                 outbasename: str, groups: dict = None):
    """
    From a dictionary of significant events, writes the 
    output files in several formats (e.g. bed files, to ggsashimi).
    
    :param defaultdict events: Events that are differentially
    spliced across all conditions
    :param list header: Header info to the main output table
    :param str outbasename: Output basename
    :param dict groups: If there are multiple comparison, this
    dict lists their names and corresponding files 
    """
    
    groups_events_map = defaultdict(list)
    strand = {"-": "minus", "+": "plus"}
    
    main_output = open(outbasename + "_rmats.tsv", 'w')
    main_output.write('\t'.join(header) + "\n")
    
    ggsashimi = open(outbasename + "_rmats_toGGsashimi.tsv", 'w')
    bed_spanning = open(outbasename + "_spanning_rmats.bed", "w")
    bed_target = open(outbasename + "_target_rmats.bed", "w")

    for coordinates_id, data in events.items():
        for i, v in enumerate(data):

            if groups is not None:
                groups_events_map[coordinates_id].append(v[-1])
                groups_with_event = ';'.join(v[-1] for v in data)

            else:
                groups_with_event = outbasename

            # If first group having this event
            if i == 0:
                main_output.write(coordinates_id + "\t" + '\t'.join(v) + "\n")
                ev_sash = [v[9], v[0], v[2], strand[coordinates_id.split(":")[1]], groups_with_event, coordinates_id]
                ggsashimi.write('\t'.join(ev_sash) + '\n')

                coord_span = data[0][-2] if groups else data[0][-1]
                _l = re.split(':|-', coord_span)
                ev_bed_span = [_l[0], _l[1], _l[2], data[0][0] + "_" + data[0][2], groups_with_event, coordinates_id]
                bed_spanning.write('\t'.join(ev_bed_span) + '\n')

                coord_target = data[0][-3] if groups else data[0][-2]
                _l = re.split(':|-', coord_target)
                ev_bed_target = [_l[0], _l[1], _l[2], data[0][0] + "_" + data[0][2], groups_with_event, coordinates_id]
                bed_target.write('\t'.join(ev_bed_target) + '\n')

            else:
                main_output.write("\t" + '\t'.join(v) + "\n")

    main_output.close()
    ggsashimi.close()
    bed_target.close()
    bed_spanning.close()

    if groups is not None:
        with open(outbasename + "_events_group_maps.csv", "w") as out:
            for k, v in groups_events_map.items():
                out.write(k + "\t" + ";".join(sorted(v)) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Script to produce readable rmats significant events between two(or more)'
                    'conditions. It picks from already filtered events through maser so'
                    'it\'s possible that by changing thresholds you wont see any'
                    'differences.')

    parser.add_argument("--input_dir", help="Path where the output rmats files are located. Use this "
                                            "argument if there is only one rmats run to be processed")
    parser.add_argument("--input_groups", help='Tab delimited file listing directories where output '
                                               'rmats files are located and corresponding comparison. '
                                               'Use this argument if multiple rmats runs need to be compared. '
                                               'Example:\ndir_to_rmats_run1 ctrl_vs_t10\n'
                                               'dir_to_rmats_run2    ctrl_vs_t20')
    parser.add_argument("-o", "--outbasename", required=True, help='Basename to the output file')
    parser.add_argument("-t", "--threshold", metavar="", type=float, default=0.2, help='dPSI threshold. Default:0.2')
    parser.add_argument("-f", "--fdr", metavar="", type=float, default=0.05, help='FDR threshold. Only events with '
                                                                                  'lower FDR value will be kept. '
                                                                                  'Default: 0.05')

    args = parser.parse_args()

    if args.input_dir and args.input_groups:
        raise ValueError("--input_dir not allowed with --input_groups. Please set one of the two options.")

    if not args.input_dir and not args.input_groups:
        raise ValueError("Please provide input data for the script using either the "
                         "--input_dir or the --input_groups argument")

    if args.input_dir:

        if os.path.isdir(args.input_dir):
            groups = None
            sign_events, header = process_rmats_files(args.input_dir,
                                                      args.threshold,
                                                      args.fdr)
        else:
            raise ValueError('Please set a valid directory where rmats file are located.')

    else:
        groups = read_groups(args.input_groups)
        sign_events, header = process_rmats_files(groups,
                                                  args.threshold,
                                                  args.fdr)

    write_output(sign_events, header, args.outbasename, groups)


if __name__ == "__main__":
    main()
