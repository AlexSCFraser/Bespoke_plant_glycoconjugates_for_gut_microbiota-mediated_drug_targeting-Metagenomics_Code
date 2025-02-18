import hashlib
import json
import math
import os
import pathlib
import pickle
import re
from codecs import strict_errors
from collections import defaultdict
from copy import deepcopy
from csv import DictReader
from csv import reader
from functools import reduce
from logging import Logger, getLogger
from inspect import getsourcefile

import matplotlib
import matplotlib.pyplot as plt

import numpy as np
import numpy
import pandas as pd
import seaborn
import seaborn as sns
import pandas

from GraphText import adjusted_labels, adjusted_names
from PULHeatmapGraphs import render_coverage_heatmap_xyloglucan, render_presence_heatmap_xyloglucan

project_folder = pathlib.Path(getsourcefile(lambda : 0)).parents[1]

ibd_metadata_files = [

]

ibd_study_name_alias_dict = {
    "HMP2_metadata": "HMP2",
}

ibd_study_id_dict = {
                 "HMP2": "HMP2"
                 }

class CazymeMetadata:
    loci_metadata: pd.DataFrame = None

    def __init__(self, pul_dataframe_file):
        self.loci_metadata = load_data_frame(pul_dataframe_file)

    def locus_ids_from_cazyme_tag(self, gene_tag: str) -> list[str]:
        try:
            locus_df = self.loci_metadata[["Gene ID", "af_note"]].dropna()\
                                                                 .query(f'af_note.str.contains("{gene_tag}")')
        except KeyError:
            locus_df = self.loci_metadata[["Gene ID", "Modularity"]].dropna()\
                                                                    .query(f'Modularity.str.contains("{gene_tag}")')

        return list(locus_df["Gene ID"].astype("string"))

    def locus_ids_from_locus_name(self, locus_name, gene_tag=None):
        try:
            locus_df = self.loci_metadata[["Gene ID", "xygul_name", "af_note"]].dropna(subset=["Gene ID", "pul_name"])\
                                                                    .query(f'xygul_name.str.startswith("{locus_name}")')
        except KeyError:
            locus_df = self.loci_metadata[["Gene ID", "pul_name", "Modularity"]].dropna(subset=["Gene ID", "pul_name"])\
                                                                    .query(f'pul_name.str.startswith("{locus_name}")')
        if gene_tag:
            if "af_note" in self.loci_metadata:
                locus_df = locus_df.dropna().query(f'af_note.str.startswith("{gene_tag}")')
            elif "Modularity" in self.loci_metadata:
                locus_df = locus_df.dropna().query(f'Modularity.str.startswith("{gene_tag}")')
            else:
                raise Exception("Locus/PUL DataFrame does not contain valid column name for Modularity!")

        return list(locus_df["Gene ID"].astype("string"))


def compute_loci_presence(run_id, loci_name, depth_threshold=1, coverage_threshold=0.7, data_folder=None) -> bool:
    if data_folder:
        filepath = os.path.join(data_folder, f"{run_id}-{loci_name}.out")
    else:
        filepath = os.path.join(os.getcwd(), "output", "97perc_test", f"{run_id}-{loci_name}.out")
    try:
        percent_below_thresh = 0.0
        with open(filepath) as csvfile:
            csv_reader = reader(csvfile, delimiter="\t")
            for line in csv_reader:
                if line[0] == "genome" and int(line[1]) < depth_threshold:
                    percent_below_thresh += float(line[4])
    except IOError:
        raise UserWarning(f"Can't read file {filepath}.")
    is_present = percent_below_thresh < (1 - coverage_threshold)
    # return percent_below_thresh < (1 - coverage_threshold)
    is_present_new = compute_loci_presence(run_id, loci_name, depth_threshold, data_folder) >= coverage_threshold #todo: WHY IS THIS LINE HERE? WHICH METHOD OF claculating is_present is correct?
    # todo: i might be in the process of refactoring this code to reduce duplication, but did not finish making above line a call to compute_loci_coverage
    assert(is_present_new == is_present)
    return is_present


def compute_loci_coverage(run_id, loci_name, locus_id_computer: CazymeMetadata, failed_list, depth_threshold=1,
                          data_folder=None, cazyme_tag=None, alternate_loci_name=None, use_loci_name_as_id=False,
                          logger: Logger = getLogger()) -> float:
    filepath_single_locus = os.path.join(data_folder, f"{run_id}-{loci_name}.out".replace(' ', '_'))
    if data_folder and not os.path.isfile(filepath_single_locus) and alternate_loci_name:
        using_alternate = True
        filepath = os.path.join(data_folder, f"{run_id}-{alternate_loci_name}.out")
    elif data_folder:
        using_alternate = False
        filepath = filepath_single_locus
    else:
        using_alternate = False
        filepath = os.path.join(os.getcwd(), "output", "97perc_test", f"{run_id}-{loci_name}.out")

    if cazyme_tag is None or cazyme_tag == "genome":
        locus_ids_old = ["genome"]
        locus_ids = locus_id_computer.locus_ids_from_locus_name(loci_name)
    else:
        locus_ids_old = locus_id_computer.locus_ids_from_cazyme_tag(cazyme_tag)
        locus_ids = locus_id_computer.locus_ids_from_locus_name(loci_name, cazyme_tag)

    if len(locus_ids) < 1:
        logger.warning(f"Locus group {loci_name} does not have any genes with tag {cazyme_tag}!")
        return math.nan

    if use_loci_name_as_id:
        gene_bases_covered_by_cazyme = defaultdict(int)
        gene_length_by_cazyme = defaultdict(int)
    else:
        gene_bases_covered_by_cazyme = dict.fromkeys(locus_ids, 0)
        gene_length_by_cazyme = dict.fromkeys(locus_ids, 0)
    try:
        if use_loci_name_as_id:
            percent_below_depth_by_cazyme = defaultdict(float)
        else:
            percent_below_depth_by_cazyme = dict.fromkeys(locus_ids, 0.0)
        percent_below_depth_old = 0.0
        with open(filepath) as csvfile:
            csv_reader = reader(csvfile, delimiter="\t")
            for line in csv_reader:
                if use_loci_name_as_id:
                    if line[0] == loci_name and int(line[1]) < depth_threshold:
                        percent_below_depth_old += float(line[4])
                    if line[0] == loci_name and int(line[1]) < depth_threshold:
                        percent_below_depth_by_cazyme[line[0]] += float(line[4])
                        gene_bases_covered_by_cazyme[line[0]] += int(line[2])
                    if line[0] == loci_name and gene_length_by_cazyme[line[0]] == 0:
                        gene_length_by_cazyme[line[0]] = int(line[3])
                else:
                    if line[0] in locus_ids_old and int(line[1]) < depth_threshold:
                        percent_below_depth_old += float(line[4])
                    if line[0] in locus_ids and int(line[1]) < depth_threshold:
                        percent_below_depth_by_cazyme[line[0]] += float(line[4])
                        gene_bases_covered_by_cazyme[line[0]] += int(line[2])
                    if line[0] in locus_ids and gene_length_by_cazyme[line[0]] == 0:
                        gene_length_by_cazyme[line[0]] = int(line[3])
    except IOError as err:
        failed_list.append(filepath)
        print(f"WARNING: Data not loaded from .out file for id: {run_id}")
        return math.nan
        # raise UserWarning(f"Can't read file {filepath}.")
    try:
        percent_below_depth = sum(gene_bases_covered_by_cazyme.values()) / sum(gene_length_by_cazyme.values())
    except ZeroDivisionError as err:
        failed_list.append(filepath)
        print(err)
        return math.nan
    if len(locus_ids) == 1 and not using_alternate:
        assert(math.isclose(percent_below_depth_old, percent_below_depth, rel_tol=1e-5))
    return 1.0 - percent_below_depth


def load_actual_read_count(data_folder, run_id, loci_name, failed_list, mixed_loci_name=None) -> int | None:
    filepath = os.path.join(data_folder, f"{run_id}-{loci_name}.out.stats")
    filepath_single_locus = os.path.join(data_folder, f"{run_id}-{loci_name}.out.stats")
    if not os.path.isfile(filepath_single_locus) and mixed_loci_name:
        using_alternate = True
        filepath = os.path.join(data_folder, f"{run_id}-{mixed_loci_name}.out.stats")
    else:
        using_alternate = False
        filepath = filepath_single_locus
    try:
        with open(filepath) as textfile:
            for line in textfile.readlines():
                if line.__contains__("paired in sequencing"):
                    paired_reads = int(line.split(' ')[0])
                elif line.__contains__("singletons"):
                    singletons = int(line.split(' ')[0])
        return int(paired_reads/2) + singletons
    except FileNotFoundError:
        failed_list.append(filepath)
        return None


def load_actual_base_count(data_folder, run_id, loci_name, failed_list, mixed_loci_name=None) -> int | None:
    filepath = os.path.join(data_folder, f"{run_id}-{loci_name}.out.basestats")
    filepath_single_locus = os.path.join(data_folder, f"{run_id}-{loci_name}.out.basestats")
    if not os.path.isfile(filepath_single_locus) and mixed_loci_name:
        using_alternate = True
        filepath = os.path.join(data_folder, f"{run_id}-{mixed_loci_name}.out.basestats")
    else:
        using_alternate = False
        filepath = filepath_single_locus
    try:
        with open(filepath) as textfile:
            for line in textfile.readlines():
                if line.__contains__("total length"):
                    bases = int(line.split('\t')[2])
        return bases
    except FileNotFoundError:
        failed_list.append(filepath)
        return math.nan


def load_metadata_file(filename):
    if not os.path.exists(filename):
        filename = os.path.join(os.getcwd(), "metadata", filename)
    data = []
    with open(filename, 'r') as csvfile:
        csv_reader = DictReader(csvfile, delimiter=',')
        # column_names = csv_reader.fieldnames
        for line in csv_reader:
            data.append(line)
    return data


def load_data_frame(filename):
    if filename.endswith("csv"):
        delimiter = ','
    elif filename.endswith("tsv"):
        delimiter = '\t'
    else:
        delimiter = ','

    return pd.read_csv(filename, delimiter=delimiter)


def count_metadata_keys(metadata):
    count_dict = defaultdict(int)
    for run in [study[0] for study in metadata.values()]:
        for key in run.keys():
            count_dict[key] += 1
    return count_dict


def get_unique_slice_indices(data_array):
    index_dict = {}

    for idx, value in enumerate(data_array):
        if type(value) == float and math.isnan(value):
            continue
        if value not in index_dict:
            index_dict[value] = (idx, idx)
        else:
            index_dict[value] = (index_dict[value][0], idx)

    return index_dict


def compute_xyloglucan_data(metadata_files, pul_list: list, data_folder=None, threshold=0.7):
    metadata = {os.path.basename(filename).split('.')[0]: load_metadata_file(filename) for filename in metadata_files}

    # sample_count = sum([len(item) for item in metadata.values()])
    sample_sizes = {}
    patient_data_by_study = {}
    study_index = 0
    sample_index = 0
    files_read = 0
    files_failed = 0
    for study_id, study_data in metadata.items():
        # study_sample_size = 0
        study_patient_ids = defaultdict(int)
        study_patient_puls = {pul_name: defaultdict(int) for pul_name in pul_list}
        study_patient_pul_coverage = {pul_name: defaultdict(list) for pul_name in pul_list}

        study_patient_theta_pul = defaultdict(int)
        study_patient_uniformis_pul = defaultdict(int)
        study_patient_fluxus_pul = defaultdict(int)
        # study_patient_pul = defaultdict(int)

        for sample_data in study_data:
            if sample_data['BioProject'] == "PRJEB6997":
                patient_identifier = sample_data["individual_name"]
            elif sample_data["BioProject"] == "PRJDB3601":
                patient_identifier = sample_data["Individual"]
            else:
                patient_identifier = sample_data["Sample Name"]
                # patient_identifier = sample_data["Run"]
            study_patient_ids[patient_identifier] += 1
            # total_pul_count = 0

            # sra_id =

            for pul_name in pul_list:
                try:
                    if compute_loci_presence(sample_data['Run'], pul_name, coverage_threshold=threshold, data_folder=data_folder): #todo: increase threshold amount
                        study_patient_puls[pul_name][patient_identifier] += 1
                    files_read += 1
                except UserWarning as error:
                    print(error.args[0])
                    files_failed += 1
                    study_patient_puls[pul_name][patient_identifier] = -1 #todo: uncomment

            # graph_data = np.ins((graph_data, np.column_stack(column)))
            sample_index += 1

        sample_sizes[study_id] = len(study_patient_ids)
        study_presence_data = {}
        for patient in study_patient_ids.keys():
            study_presence_data[patient] = [study_index
                                            # min(study_patient_fluxus_pul[patient], 1),
                                            # min(study_patient_theta_pul[patient], 1),
                                            # min(study_patient_uniformis_pul[patient], 1)
                                            ]
            for pul_name in pul_list:
                # if patient in study_patient_puls[pul_name]:
                study_presence_data[patient].append(min(study_patient_puls[pul_name][patient], 1))
                # else:
                #     study_presence_data[patient].append(-1)
            study_presence_data[patient].append(sum(study_presence_data[patient][1:]))
        patient_data_by_study[study_id] = study_presence_data
        print("Study analyzed")
        study_index += 1
    return patient_data_by_study, metadata, files_failed

def compute_xyloglucan_data_df(metadata_files, computed_df_folder, pul_list: list, pul_dataframe_file, pul_id: str,
                               data_folder=None, coverage_threshold=0.7, cazyme_list: tuple[str, ...] = "genome",
                               pul_alias_dict: dict = {}, keep_patient_timeseries=False, mixed_loci_name=None,
                               use_loci_name_as_id=False):
    opts = ''.join(metadata_files) + ''.join(cazyme_list) + str(keep_patient_timeseries) + str(coverage_threshold) + \
           ''.join(pul_list) + data_folder + str(use_loci_name_as_id) + pul_id
    computed_df_file = hashlib.md5(opts.encode('utf-8')).hexdigest() + '.pickle'
    picklepath = os.path.join(computed_df_folder, computed_df_file)
    if os.path.isfile(picklepath):
        with open(picklepath, 'rb') as file:
            loaded_data = pickle.load(file)
        return loaded_data["metadata"], loaded_data["files_failed"], loaded_data["pul_list"], loaded_data["coverage_fieldnames"], \
            loaded_data["threshold_fieldnames"]

    metadata: dict[str: pd.DataFrame] = {pathlib.Path(filename).stem: load_data_frame(filename) for filename in metadata_files}

    file_failed_list = []
    for study_id, study_data in metadata.items():
        coverage_fieldnames = defaultdict(list)
        threshold_fieldnames = defaultdict(list)
        required_fieldname_mapping = {"patient_id": ["Individual", "Sample_Name", "sample_name", "Data_ID",
                                                     "individual_name", "Subject", "subject_accession",
                                                     "submitted_subject_id", "BioSample"],
                                      "study_name": ["BioProject", "study_id"],
                                      "sra_num": ["Run", "sample_accession"],
                                      "phenotype": ["Phenotype", "disease"],
                                      "environment_material": ["environment_(material)", "Isolation_source", "env_material", "sample_type"]
                                      }
        optional_fieldname_mapping = {"sex": ["Sex", "Gender", "gender"],
                                      "age": ["Age", "Age_x", "Age_y", "age_category"],
                                      "Bases": [],
                                      "bmi": ["BMI"],
                                      "timepoint": [],
                                      "immunosuppressants": [],
                                      "steroids": [],
                                      "antibiotics": [],
                                      "five_asa": ["5ASA"],
                                      "Organism": []
                                      }
        data_replacements = {
            "Faeces": "Stool",
            "faeces": "Stool",
            "Feces": "Stool",
            "feces": "Stool",
            "Human feces": "Stool",
            "G_DNA_Stool": "Stool",
            "Crohn's Disease": "CD",
            "Healthy": "Control",
            "HC": "Control",
            "False": False,
            "True": True,
            "None": pandas.NA,
        }

        # Clean up data
        for original, replacement in data_replacements.items():
            study_data = study_data.replace(to_replace=original, value=replacement)
        # study_data = study_data.replace(to_replace="Crohn's Disease", value="CD")
        # study_data = study_data.replace(to_replace="Healthy", value="Control")
        # study_data = study_data.replace(to_replace="HC", value="Control")
        # study_data = study_data.replace(to_replace="False", value=False)
        # study_data = study_data.replace(to_replace="True", value=True)
        # possibly change pandas.NA to numpy.NaN, since pandas.NA seems bugged
        # study_data = study_data.replace(to_replace="None", value=pandas.NA)

        if "Phenotype" in study_data:
            study_data["Phenotype"] = study_data["Phenotype"].replace(np.nan, "NA")
        study_data = study_data.dropna(axis="columns", how="all")

        columns_to_keep = list(required_fieldname_mapping.keys())
        # columns_to_keep.update(optional_fieldname_mapping.keys())
        if not keep_patient_timeseries:
            columns_to_keep.remove("patient_id")

        for target_fieldname, possible_fieldnames in required_fieldname_mapping.items():
            if target_fieldname not in study_data:
                # possible_patient_fieldnames = possible_fieldnames
                for potential_id in possible_fieldnames:
                    if potential_id in study_data:
                        study_data[target_fieldname] = study_data[potential_id]
                        break
            if target_fieldname not in study_data:
                raise Exception(f"No valid {target_fieldname} field in dataset from file: {study_id}")

        for target_fieldname, possible_fieldnames in optional_fieldname_mapping.items():
            if target_fieldname not in study_data:
                # possible_patient_fieldnames = possible_fieldnames
                for potential_id in possible_fieldnames:
                    if potential_id in study_data:
                        study_data[target_fieldname] = study_data[potential_id]
                        break
            if target_fieldname in study_data:
                columns_to_keep.append(target_fieldname)
            else:
                print(f"WARNING: No valid {target_fieldname} field in dataset from file: {study_id}")

        # load coverage data into dataframe and calculate PUL presence
        locus_id_computer = CazymeMetadata(pul_dataframe_file)
        # study_data = study_data.assign(xygul_count=lambda initial_count: 0) # replaced by simpler below
        study_data[f'{pul_id.lower()}_count'] = 0
        for cazyme in cazyme_list:
            for pul_name in pul_list:
                if cazyme != "genome" and len(locus_id_computer.locus_ids_from_locus_name(pul_name, cazyme)) < 1:
                    continue  # skips trying to calculate coverage for cazymes which aren't in specific PULs
                newline = '\n'
                formatted_pul_name = f'{pul_alias_dict[pul_name] if pul_name in pul_alias_dict else pul_name}' \
                                     f'{"" if cazyme == "genome" or cazyme is None else newline + cazyme}'
                # Calculate total locus coverage
                coverage_fieldname = formatted_pul_name + "\ncoverage"
                columns_to_keep.append(coverage_fieldname)
                coverage_fieldnames[cazyme].append(coverage_fieldname)
                study_data[coverage_fieldname] = study_data['sra_num'].apply(lambda sra_num: compute_loci_coverage(
                    sra_num, pul_name, locus_id_computer, file_failed_list, data_folder=data_folder,
                    cazyme_tag=cazyme, alternate_loci_name=mixed_loci_name,
                    use_loci_name_as_id=use_loci_name_as_id))
                # calculate whether the current PUL locus has coverage exceeding threshold and store result in new
                # column
                threshold_fieldname = formatted_pul_name + "\npresent"
                columns_to_keep.append(threshold_fieldname)
                threshold_fieldnames[cazyme].append(threshold_fieldname)
                study_data[threshold_fieldname] = study_data[coverage_fieldname].apply(
                                                    lambda coverage: 1 if coverage >= coverage_threshold else 0)

        # load actual read counts
        # study_data["actual_reads"] = study_data['sra_num']\
        #     .apply(lambda sra_num: load_actual_read_count(data_folder, sra_num, pul_name, file_failed_list,
        #                                                   mixed_loci_name))
        study_data["actual_bases"] = study_data['sra_num']\
            .apply(lambda sra_num: load_actual_base_count(data_folder, sra_num, pul_name, file_failed_list,
                                                          mixed_loci_name))

        columns_to_keep.append("actual_bases")
        # todo: reimplement assertion after loading actual base count from basestats file
        if "Bases" in study_data:
            # todo reimplement assertion later once twin study data issue is solved
            # assert((study_data["actual_bases"] == study_data["Bases"]).all())
            pass
        else:
            study_data["Bases"] = study_data["actual_bases"]
            columns_to_keep.append("Bases")
        # assert((study_data["actual_read_base_estimate"] / study_data["Bases"] < 0.01).all())

        max_cols = reduce(lambda a, b: a+b, coverage_fieldnames.values()) + \
                   reduce(lambda a, b: a+b, threshold_fieldnames.values()) + \
                   ["Bases", "actual_bases"]
        if keep_patient_timeseries:
            # handle duplicate patient+timepoint combination by selecting the highest coverage of each PUL
            # each timepoint per patient
            # todo: FIX and test new version of timeseries data loading works
            max_df = study_data[["patient_id", "timepoint"] + max_cols].groupby(["patient_id", "timepoint"]).max()
            study_data = max_df.join(study_data[list(set(study_data.columns) - set(max_cols))].set_index(["patient_id", "timepoint"])).reset_index()
            # old version
            # study_data = study_data.groupby(["patient_id", "timepoint"], as_index=False).max()[columns_to_keep]
        elif not keep_patient_timeseries:
            # handle duplicate patients by selecting the highest coverage of each PUL each observation
            max_df = study_data[["patient_id"] + max_cols].groupby("patient_id").max()
            sample_cols = ["actual_bases", "Bases", "sample_id", "timepoint", "sra_num", "sample_accession"]
            cols_to_drop = list(filter(lambda a: a in study_data, max_cols + sample_cols))
            relevant_data = study_data.drop(cols_to_drop, axis=1) \
                                      .drop_duplicates(subset=['patient_id'], keep='first') \
                                      .dropna(subset=['patient_id']) \
                                      .set_index("patient_id")
            assert(len(relevant_data) == len(max_df))
            study_data = max_df.join(relevant_data).reset_index()

        # sum PUL counts and add presence column
        for cazyme in cazyme_list:
            count_attribute = f"{pul_id.lower()}{'' if cazyme == 'genome' else '_' + cazyme}_count"
            any_pul_attribute = f"Any\n{pul_id if cazyme == 'genome' else cazyme}\npresent"
            study_data[count_attribute] = study_data[threshold_fieldnames[cazyme]].agg('sum', axis=1)
            study_data[any_pul_attribute] = study_data[count_attribute].apply(lambda num: 1 if num >= 1 else 0)
            threshold_fieldnames[cazyme].append(any_pul_attribute)

            max_coverage_fieldname = f"max_{pul_id.lower() if cazyme == 'genome' or cazyme is None else cazyme}_coverage"
            study_data[max_coverage_fieldname] = study_data[coverage_fieldnames[cazyme]].max(1)
            coverage_fieldnames[cazyme].insert(0, max_coverage_fieldname)

        # add column to indicate whether treatments were received
        treatments = ["antibiotics", "five_asa", "steroids", "immunosuppressants"]
        treatments_in_dataset = list(set.intersection(set(treatments), set(study_data.keys())))
        if treatments_in_dataset:
            study_data["any_treatment"] = study_data[treatments_in_dataset]\
                .replace(to_replace=numpy.NaN, value=False)\
                .astype(bool)\
                .agg("max", axis="columns")

        if "antibiotics" in study_data and "steroids" in study_data:
            study_data["antibiotics_or_steroids"] = study_data[["antibiotics", "steroids"]]\
                .replace(to_replace="None", value=pandas.NA)\
                .agg("max", axis="columns")

        # add base pairs in gigabases
        study_data["GigaBases"] = study_data["Bases"] / 1e9

        metadata[study_id] = study_data

    file_failed_list = list(set(file_failed_list))
    files_failed = len(file_failed_list)
    if files_failed > 0:
        print(f"List of failed files: ", '\n\t'.join(file_failed_list))
    else:
        pass
        data_to_serialize = {"metadata": metadata,
                             "files_failed": files_failed,
                             "pul_list": pul_list,
                             "coverage_fieldnames": coverage_fieldnames,
                             "threshold_fieldnames": threshold_fieldnames}
        os.makedirs(os.path.dirname(picklepath), exist_ok=True)
        with open(picklepath, 'wb') as file:
            pickle.dump(data_to_serialize, file)

    failed_id_list = list(set([pathlib.PurePath(fn).stem.split('-')[0] for fn in file_failed_list]))
    if len(failed_id_list) > 0:
        print(failed_id_list)

    return metadata, files_failed, pul_list, coverage_fieldnames, threshold_fieldnames


def render_heatmap(patient_data_by_study, metadata, study_xlabels):
    total_patient_count = 0
    for study in patient_data_by_study.values():
        total_patient_count += len(study)
    graph_data = numpy.full((5, total_patient_count), -1, dtype=numpy.int8)
    total_patient_index = 0
    fluxus_count = 0
    theta_count = 0
    uniformis_count = 0
    any_count = 0
    for study_id, study_data in patient_data_by_study.items():
        for patient in study_data:
            graph_data[0, total_patient_index] = patient_data_by_study[study_id][patient][0]
            graph_data[1, total_patient_index] = patient_data_by_study[study_id][patient][1]
            graph_data[2, total_patient_index] = patient_data_by_study[study_id][patient][2]
            graph_data[3, total_patient_index] = patient_data_by_study[study_id][patient][3]
            graph_data[4, total_patient_index] = patient_data_by_study[study_id][patient][4]
            if graph_data[1, total_patient_index] > 0:
                fluxus_count += 1
            if graph_data[2, total_patient_index] > 0:
                theta_count += 1
            if graph_data[3, total_patient_index] > 0:
                uniformis_count += 1
            if graph_data[4, total_patient_index] > 0:
                any_count += 1
            total_patient_index += 1

    fluxus_freq = fluxus_count / total_patient_index
    theta_freq = theta_count / total_patient_index
    uniformis_freq = uniformis_count / total_patient_index
    any_freq = any_count / total_patient_index

    figure, axis = matplotlib.pyplot.subplots(2, 1, sharex=False, gridspec_kw={'height_ratios': [1, 50]})
    # figure.set_figure_size()
    # matplotlib.pyplot
    # axis[0].set_autoscaley_on()
    ax = seaborn.heatmap(graph_data[0:1, :], ax=axis[0], vmin=0, vmax=11,
                         cmap=["red", "yellow", "green",
                               "cyan", "blue", "violet", "red", "green", "orangered", "orange", "salmon",
                               "green"], cbar=False
                         )
    study_label_pos = []
    curr_pos = 0
    for study in patient_data_by_study.values():
        study_label_pos.append(curr_pos + len(study)/2)
        curr_pos += len(study)
    ax.set_xticks(study_label_pos)
    ax.set_xticklabels(study_xlabels, size=5)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, pad=10, labelrotation=45)
    # matplotlib.axes.Axes.set_aspect(axis[1], 1)
    ax2 = seaborn.heatmap(graph_data[1:5, :], ax=axis[1], vmin=-1, vmax=3,
                          cmap=["purple", "black", "blue", "yellow", "red"],
                          cbar_kws={"location": "bottom", "fraction": 0.05, "label": "PUL count",
                                    "ticks": [0.25, 1, 1.8, 2.6], "format": "%1.0f"},
                          yticklabels=["B. fluxus 1,3GUL", "B. thetaiotaomicron 1,3GUL", "B. uniformis 1,3GUL",
                                       "All species 1,3GUL"], xticklabels=[]
                          )
    ax2.set_yticklabels(["B. fluxus\n1,3GUL", "B. thetaiotaomicron\n1,3GUL", "B. uniformis\n1,3GUL",
                         "All species\n1,3GUL"], size=5)
    percent_ax = ax2.twinx()
    percent_ax.set_yticklabels([f"{fluxus_freq*100: 0.2f}", f"{theta_freq*100: 0.2f}", f"{uniformis_freq*100: 0.2f}",
                                f"{any_freq*100: 0.2f}"], size=6)
    percent_ax.set_yticks(numpy.arange(4)+0.5)
    percent_ax.set_ylim(ax2.get_ylim())
    percent_ax.set_ylabel("Frequencies (%)")

    # ax2.set_aspect(1)
    # aspect = 1
    # im = ax2.get_images()
    # extent = im[0].get_extent()
    # ax2.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    ax2.figure.set_dpi(300)
    # ax = seaborn.heatmap(graph_data, cmap=["black", "blue", "yellow", "red", "red", "yellow", "green",
    #                                        "teal", "blue", "purple", "red", "green", "orangered", "orange", "green"])
    # giving title to the plot
    matplotlib.pyplot.title('Heatmap')

    # function to show plot
    matplotlib.pyplot.show()
    print("yo")


def render_heatmap_xyloglucan(patient_data_by_study, metadata, study_xlabels, pul_list, pul_names, title):
    total_patient_count = 0
    for study in patient_data_by_study.values():
        total_patient_count += len(study)
    # numby array contains rows for each PUL, plus the first row for study id, plus the last row for total pul count
    graph_data = numpy.full((len(pul_list) + 2, total_patient_count), -1, dtype=numpy.int8)
    total_patient_index = 0
    species_counts = [0] * (len(pul_list) + 1)  # add one to length to hold counts/freqs for any species
    species_freqs = [0] * (len(pul_list) + 1)
    # fluxus_count = 0
    # theta_count = 0
    # uniformis_count = 0
    any_count = 0
    for study_id, study_data in patient_data_by_study.items():
        for patient in study_data:
            # graph_data[0, total_patient_index] = patient_data_by_study[study_id][patient][0]
            # graph_data[1, total_patient_index] = patient_data_by_study[study_id][patient][1]
            # graph_data[2, total_patient_index] = patient_data_by_study[study_id][patient][2]
            # graph_data[3, total_patient_index] = patient_data_by_study[study_id][patient][3]
            # graph_data[4, total_patient_index] = patient_data_by_study[study_id][patient][4]
            graph_data[:, total_patient_index] = patient_data_by_study[study_id][patient]

            for i in range(len(pul_list)+1):
                if graph_data[i+1, total_patient_index] > 0:
                    species_counts[i] += 1

            # if graph_data[1, total_patient_index] > 0:
            #     fluxus_count += 1
            # if graph_data[2, total_patient_index] > 0:
            #     theta_count += 1
            # if graph_data[3, total_patient_index] > 0:
            #     uniformis_count += 1
            # if graph_data[-1, total_patient_index] > 0:
            #     any_count += 1

            total_patient_index += 1

    species_freqs = list(map(lambda x: f"{x/total_patient_count: 0.2f}", species_counts))
    # fluxus_freq = fluxus_count / total_patient_index
    # theta_freq = theta_count / total_patient_index
    # uniformis_freq = uniformis_count / total_patient_index
    # any_freq = any_count / total_patient_index

    figure, axis = matplotlib.pyplot.subplots(2, 1, sharex=False, gridspec_kw={'height_ratios': [1, 50]})
    # figure.set_figure_size()
    # matplotlib.pyplot
    # axis[0].set_autoscaley_on()
    ax = seaborn.heatmap(graph_data[0:1, :], ax=axis[0], vmin=0, vmax=11,
                         cmap=["red", "yellow", "green",
                               "cyan", "blue", "violet", "red", "green", "orangered", "orange", "salmon",
                               "green"], cbar=False
                         )
    study_label_pos = []
    curr_pos = 0
    for study in patient_data_by_study.values():
        study_label_pos.append(curr_pos + len(study)/2)
        curr_pos += len(study)
    ax.set_xticks(study_label_pos)
    ax.set_xticklabels(study_xlabels, size=5)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, pad=10, labelrotation=25)
    # matplotlib.axes.Axes.set_aspect(axis[1], 1)
    ax2 = seaborn.heatmap(graph_data[1:, :], ax=axis[1], vmin=-1, vmax=4,
                          cmap=["purple", "black", "blue", "green", "yellow", "red"],
                          cbar_kws={"location": "bottom", "fraction": 0.05, "label": "PUL count",
                                    "ticks": [0.25, 1, 1.8, 2.6, 3.6], "format": "%1.0f"},
                          yticklabels=pul_names + ["Species Count"], xticklabels=[]
                          )
    ax2.set_yticklabels(pul_names + ["XyGUL Count"], size=5)
    percent_ax = ax2.twinx()
    percent_ax.set_yticklabels(species_freqs, size=6)
    percent_ax.set_yticks(numpy.arange(len(pul_list)+1)+0.5)
    percent_ax.set_ylim(ax2.get_ylim())
    percent_ax.set_ylabel("Frequencies (%)")

    # ax2.set_aspect(1)
    # aspect = 1
    # im = ax2.get_images()
    # extent = im[0].get_extent()
    # ax2.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    ax2.figure.set_dpi(300)
    # ax = seaborn.heatmap(graph_data, cmap=["black", "blue", "yellow", "red", "red", "yellow", "green",
    #                                        "teal", "blue", "purple", "red", "green", "orangered", "orange", "green"])
    # giving title to the plot
    matplotlib.pyplot.title(title)

    # function to show plot
    matplotlib.pyplot.show()
    print("yo")


def merge_dataframes_with_seperator_row(dataframes: list[pd.DataFrame], num_empty_rows=1):
    empty_row = {col: [math.nan] for col in dataframes[0].columns}
    index_name = dataframes[0].index.name
    if index_name is None:
        index_names = dataframes[0].index.names
    else:
        index_names = [index_name]
    for name in index_names:
        empty_row[name] = '0'
    empty_row_df = pd.DataFrame(empty_row)
    empty_row_df = empty_row_df.set_index(index_names)
    merged = reduce(lambda df_1, df_2:
                    pd.concat([df_1] + [empty_row_df]*num_empty_rows + [df_2]), dataframes)
    return merged


def render_coverage_comparison(dataframe: pd.DataFrame, col1: str, col2: str, out_folder: os.PathLike | str,
                               title: str = None, thresh: float = 0.7):

    point_alpha = 0.35

    g = sns.jointplot(x=col1, y=col2, data=dataframe,

                      xlim=(0, 1.0), ylim=(0, 1.0),
                      color="orange", alpha=0, ratio=4)

    g.ax_marg_x.set_title(title, pad=-4)
    g.ax_joint.axvline(thresh, color="black", linestyle="--", linewidth=1)
    g.ax_joint.axhline(thresh, color="black", linestyle="--", linewidth=1)

    quadrant1 = dataframe[(dataframe[col1] >= thresh) & (dataframe[col2] >= thresh)]
    quadrant2 = dataframe[(dataframe[col1] < thresh) & (dataframe[col2] >= thresh)]
    quadrant3 = dataframe[(dataframe[col1] < thresh) & (dataframe[col2] < thresh)]
    quadrant4 = dataframe[(dataframe[col1] >= thresh) & (dataframe[col2] < thresh)]

    col1_species = col1.split('\n')[0]
    col1_leg_label = '\n'.join(col1_species.split(' ')) + ' only'
    col2_species = col2.split('\n')[0]
    col2_leg_label = '\n'.join(col2.split('\n')[0].split(' ')) + ' only'

    g.ax_joint.scatter(quadrant1[col1], quadrant1[col2], color='purple', label='Both\npresent', alpha=point_alpha)
    g.ax_joint.scatter(quadrant4[col1], quadrant4[col2], color='blue', label=col1_leg_label, alpha=point_alpha)
    g.ax_joint.scatter(quadrant2[col1], quadrant2[col2], color='red', label=col2_leg_label, alpha=point_alpha)
    g.ax_joint.scatter(quadrant3[col1], quadrant3[col2], color='green', label=f'Neither\npresent', alpha=point_alpha)

    g.ax_joint.legend(loc='upper left', bbox_to_anchor=(1, 1.29), fontsize=6)
    filename = f"{title + ' ' if title else ''}{col1_species} {col2_species} {thresh:.1f} coverage comparison"
    plt.savefig(os.path.join(out_folder, f"{filename}.png"), dpi=300)
    plt.savefig(os.path.join(out_folder, f"{filename}.eps"), dpi=300)
    plt.show()
    pass


def render_coverage_comparison_regression(dataframe: pd.DataFrame, col1: str, col2: str,
                                          out_folder: os.PathLike | str, title: str = None, thresh: float = 0.7):
    point_alpha = 0.35
    corr_coef = np.corrcoef(dataframe[col1], dataframe[col2])[0, 1]
    r_squared = corr_coef ** 2


    g = sns.jointplot(x=col1, y=col2, data=dataframe,
                      kind="reg",
                      xlim=(0, 1.0), ylim=(0, 1.0),
                      color="teal", ratio=4,
                      scatter_kws={"alpha": point_alpha},
                      line_kws=dict(color="black"),
                      ci=None)

    g.ax_joint.axvline(thresh, color="black", linestyle="--", linewidth=1)
    g.ax_joint.axhline(thresh, color="black", linestyle="--", linewidth=1)
    g.ax_joint.annotate(f'R = {corr_coef:.2f}',
                        xy=(1.04, 1.10), xycoords='axes fraction',
                        fontsize=12, ha='left', va='top')
    g.ax_joint.annotate(f'$R^2$ = {r_squared:.2f}',
                        xy=(1.04, 1.15), xycoords='axes fraction',
                        fontsize=12, ha='left', va='top')
    g.ax_marg_x.set_title(title, pad=-4)

    col1_species = col1.split('\n')[0]
    col2_species = col2.split('\n')[0]
    filename = f"{title + ' ' if title else ''}{col1_species} {col2_species} {thresh:.1f} coverage regression"
    plt.savefig(os.path.join(out_folder, f"{filename}.png"), dpi=300)
    plt.savefig(os.path.join(out_folder, f"{filename}.eps"), dpi=300)
    plt.show()
    pass


def render_stacked_bar_comparison(dataframe, title, subtitle, coverage_threshold, col_label, data_label,
                                  id_ordering, first_colname, second_colname,graphics_output_folder):
    dataframe["stack_category"] = None

    dataframe.loc[(dataframe[first_colname] < coverage_threshold) & (dataframe[second_colname] < coverage_threshold), "stack_category"] = "Neither present"
    dataframe.loc[(dataframe[first_colname] < coverage_threshold) & (dataframe[second_colname] >= coverage_threshold), "stack_category"] = second_colname + " present"
    dataframe.loc[(dataframe[first_colname] >= coverage_threshold) & (dataframe[second_colname] < coverage_threshold), "stack_category"] = first_colname + " present"
    dataframe.loc[(dataframe[first_colname] >= coverage_threshold) & (dataframe[second_colname] >= coverage_threshold), "stack_category"] = "Both present"

    presence_cat_ordering = ["Neither present",
                             second_colname + " present",
                             first_colname + " present",
                             "Both present"]

    palette = [sns.color_palette()[9], sns.color_palette()[8], sns.color_palette("hls", 8)[3], sns.color_palette("hls", 8)[5]]

    # Melt data into long format
    tmp_df = dataframe.reset_index(level="patient_id")[["stack_category"]].reset_index()
    proportions = pd.crosstab(tmp_df['study_name'], tmp_df['stack_category'], normalize='index')
    proportions = proportions.reindex(index=id_ordering, columns=presence_cat_ordering)
    df_melt = proportions.reset_index().melt(id_vars=[col_label], var_name='Presence Category', value_name='Coverage')

    # render stacked bar chart
    ax = sns.histplot(data=df_melt, x=col_label, hue='Presence Category', weights='Coverage', multiple='stack',
                      palette = palette)
    plt.tick_params(axis='x', labelsize=6)
    ax.set_ylabel("Proportion of Individuals in Category")
    ax.set_title(title)
    plt.suptitle(title)

    file_stem = f"{title} {subtitle}"
    plt.savefig(os.path.join(graphics_output_folder, f"{file_stem}.png"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{file_stem}.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{file_stem}.svg"), dpi=300)
    plt.show()

    pass


def graph_rendering_pipeline(dataframe_dict, study_alias_dict, study_id_dict, pul_list, pul_names, coverage_attributes,
                             presence_attributes, cazyme_name, coverage_threshold: float, alignment_threshold: float,
                             graphics_base_folder, pul_tag, group_by_study, vertical_group_labels, study_id_mapper=None,
                             coverage_comparisons: list[tuple] = None):
    seaborn.set(font="monospace")

    for k, v in dataframe_dict.items():
        dataframe_dict[k] = dataframe_dict[k].replace(study_id_mapper)

    treatments = ["antibiotics", "five_asa", "steroids", "immunosuppressants", "any_treatment",
                  "antibiotics_or_steroids"]
    heatmap_index_levels = ["study_name", "phenotype", "patient_id"] if group_by_study else ["phenotype", "study_name", "patient_id"]
    mapper_insert_1 = "" if cazyme_name == 'genome' or cazyme_name == '' else cazyme_name + '_'
    mapper_insert_2 = "" if cazyme_name == 'genome' or cazyme_name == '' else cazyme_name + '\n'
    mapper_insert_3 = pul_tag.lower() if cazyme_name == 'genome' or cazyme_name == '' else cazyme_name
    any_pul_present = presence_attributes[-1]
    # todo: change column remappers to be parameter or computed from parameter
    presence_column_mapper = {f"xygul_{mapper_insert_1}count": 'Count',
                              any_pul_present: "Any",
                              f'B. ovatus\n{mapper_insert_2}present': 'B. ovatus',
                              f'B.\ncellulosilyticus\n{mapper_insert_2}present': 'B. cell.',
                              f'B. uniformis\nxyGUL 1\n{mapper_insert_2}present': 'B. uniformis\nxyGUL 1',
                              f'B. uniformis\nxyGUL 2\n{mapper_insert_2}present': 'B. uniformis\nxyGUL 2',
                              f'B. fluxus\n{mapper_insert_2}present': 'B. fluxus'}
    coverage_column_mapper = {f"max_{mapper_insert_3}_coverage": 'Max',
                              f'B. ovatus\n{mapper_insert_2}coverage': 'B. ovatus',
                              f'B.\ncellulosilyticus\n{mapper_insert_2}coverage': 'B. cell.',
                              f'B. uniformis\nxyGUL 1\n{mapper_insert_2}coverage': 'B. uniformis\nxyGUL 1',
                              f'B. uniformis\nxyGUL 2\n{mapper_insert_2}coverage': 'B. uniformis\nxyGUL 2',
                              f'B. fluxus\n{mapper_insert_2}coverage': 'B. fluxus'}
    merged_dataframe = pd.concat(dataframe_dict, names=["filename"]) \
        .droplevel("filename") \
        .reset_index() \
        .replace(study_alias_dict) \
        .set_index(["study_name", "patient_id"])

    merged_dataframe_alt_grouping = merge_dataframes_with_seperator_row(list(dataframe_dict.values()), 3)
    graphics_output_folder = os.path.join(graphics_base_folder,
                                          f"coverage-{coverage_threshold}_alignment-{alignment_threshold}")
    if not os.path.isdir(graphics_output_folder):
        os.makedirs(graphics_output_folder)

    count_name = f"{pul_tag.lower()}{'' if cazyme_name == 'genome' or cazyme_name == '' else '_' + cazyme_name}_count"
    # heatmap_presence_attributes = ["patient_id", "study_id", "phenotype"]
    # heatmap_presence_attributes = ["study_id", count_name]
    heatmap_presence_attributes = [count_name]
    stratify_properties = ["phenotype"]
    # name_prefix = '' if cazyme_name == 'genome' or cazyme_name == '' else cazyme_name + '_'
    name_prefix = ''
    name_title_insert = pul_tag if cazyme_name == 'genome' or cazyme_name is None or cazyme_name == '' else cazyme_name
    align_thresh_insert = f"{alignment_threshold * 100: 2.0f}% percent identity threshold"
    coverage_thresh_insert = f"{coverage_threshold * 100: 2.0f}% coverage threshold"
    align_cov_subtitle = f"{align_thresh_insert}, {coverage_thresh_insert}"
    ordering = ["Control", "UC", "CD", "NA"]

    merged_dataframe_plot_presence = merged_dataframe[stratify_properties + heatmap_presence_attributes + presence_attributes]
    merged_dataframe_plot_presence.pop(any_pul_present)
    merged_dataframe_plot_presence = merged_dataframe_plot_presence.set_index("phenotype", append=True).reorder_levels(heatmap_index_levels)
    if not group_by_study:
        merged_dataframe_plot_presence = merged_dataframe_plot_presence.reindex(ordering, level=0)
    merged_dataframe_plot_presence.rename(columns=presence_column_mapper, inplace=True)
    # line_nums = heatmap_data[max(heatmap_data.columns[0:-1])].groupby(level=0).count()
    # todo: uncomment below
    # render_bases_coverage_scatterplot(merged_dataframe, cazyme_name, graphics_output_folder, name_prefix,
    #                                   title=f"{name_title_insert} coverage by sequencing depth, all studies",
    #                                   subtitle=align_thresh_insert + ", logistic regression", ordering=ordering)
    #
    # merged_dataframe_antibiotics_removed = merged_dataframe[merged_dataframe.antibiotics != True]
    # render_bases_coverage_scatterplot(merged_dataframe_antibiotics_removed, cazyme_name, graphics_output_folder,
    #                                   "Antibiotics removed, " + name_prefix,
    #                                   title=f"{name_title_insert} coverage by sequencing depth, all studies",
    #                                   subtitle=align_thresh_insert + ", logistic regression", ordering=ordering)
    #
    # csprism_dataframe_antibiotics_removed = merged_dataframe_antibiotics_removed[merged_dataframe_antibiotics_removed.study_name == "CS-PRISM"]
    # render_bases_coverage_scatterplot(csprism_dataframe_antibiotics_removed, cazyme_name, graphics_output_folder,
    #                                   "CS-PRISM, antibiotics removed, " + name_prefix,
    #                                   title=f"{name_title_insert} coverage by sequencing depth, CS-PRISM",
    #                                   subtitle=align_thresh_insert + ", logistic regression", ordering=ordering)
    #
    # csprism_dataframe_any_removed = merged_dataframe.query("any_treatment == False and study_name == 'CS-PRISM'")
    # render_bases_coverage_scatterplot(csprism_dataframe_any_removed, cazyme_name, graphics_output_folder,
    #                                   "CS-PRISM, no treatments, " + name_prefix,
    #                                   title=f"{name_title_insert} coverage by sequencing depth, CS-PRISM",
    #                                   subtitle=align_thresh_insert + ", logistic regression", ordering=ordering)
    #
    # hmp2_dataframe_antibiotics_removed = merged_dataframe_antibiotics_removed[merged_dataframe_antibiotics_removed.study_name == "HMP2"]
    # render_bases_coverage_scatterplot(hmp2_dataframe_antibiotics_removed, cazyme_name, graphics_output_folder,
    #                                   "HMP2, antibiotics removed, " + name_prefix,
    #                                   title=f"{name_title_insert} coverage by sequencing depth, HMP2",
    #                                   subtitle=align_thresh_insert + ", logistic regression", ordering=ordering)

    render_bases_coverage_by_study(merged_dataframe.reset_index(level="study_name"), cazyme_name, graphics_output_folder, name_prefix, study_id_dict,
                                   title="", subtitle="")

    render_presence_heatmap_xyloglucan(merged_dataframe_plot_presence, align_cov_subtitle,
                                       name_title_insert + " Heatmap", alignment_threshold, coverage_threshold,
                                       name_prefix, graphics_output_folder)

    merged_dataframe_plot_coverage = merged_dataframe[coverage_attributes]
    if not group_by_study:
        merged_dataframe_plot_coverage = merged_dataframe[stratify_properties + coverage_attributes]
        merged_dataframe_plot_coverage = merged_dataframe_plot_coverage.set_index("phenotype", append=True) \
            .reorder_levels(heatmap_index_levels) \
            .reindex(ordering, level=0)

    merged_dataframe_plot_coverage.rename(columns=coverage_column_mapper, inplace=True)

    merged_dataframe_alt_grouping["dummy_var"] = "dummy"
    merged_dataframe_alt_grouping = merged_dataframe_alt_grouping.reset_index() \
        .set_index(['dummy_var', 'study_name', 'patient_id'])
    fixed_column_names = adjusted_names(list(merged_dataframe_alt_grouping[coverage_attributes].columns),
                                        name_title_insert=name_title_insert)
    fixed_column_names = list(map(lambda string: string.replace(' ', '\n'), fixed_column_names))
    y_ticks_labels = {k: np.mean(v) for k, v in get_unique_slice_indices(merged_dataframe_alt_grouping.index.get_level_values(1).to_numpy()).items()}
    y_ticks = np.array(list(y_ticks_labels.values()))
    render_coverage_heatmap_xyloglucan(alignment_threshold, coverage_attributes, graphics_output_folder,
                                       merged_dataframe_alt_grouping[coverage_attributes],
                                       name_prefix, name_title_insert,
                                       fixed_column_names, align_thresh_insert,
                                       f"{name_title_insert} Coverage Heatmap", alignment_threshold,
                                       vertical_group_labels=vertical_group_labels,
                                       yticks=y_ticks, ytick_labels=y_ticks_labels.keys())

    # todo: temp for MLG paper, remove later
    reordered_coverage_attributes = coverage_attributes[1:]
    reordered_fixed_column_names = fixed_column_names[1:]
    try:
        reordered_coverage_attributes.remove("Segatella buccae\ncoverage")
        reordered_coverage_attributes.remove('Hoylesella loescheii\ncoverage')
        reordered_coverage_attributes.remove('Bacteroides ovatus\ncoverage')
        reordered_coverage_attributes.append('Bacteroides ovatus\ncoverage')

        reordered_fixed_column_names.remove("Segatella\nbuccae")
        reordered_fixed_column_names.remove('Hoylesella\nloescheii')
        reordered_fixed_column_names.remove('Bacteroides\novatus')
        reordered_fixed_column_names.append('Bacteroides\novatus')
    except:
        pass

    render_coverage_heatmap_xyloglucan(alignment_threshold, reordered_coverage_attributes, graphics_output_folder,
                                       merged_dataframe_alt_grouping[reordered_coverage_attributes],
                                       name_prefix, name_title_insert,
                                       reordered_fixed_column_names, align_thresh_insert,
                                       f"{name_title_insert} Coverage Heatmap", alignment_threshold,
                                       vertical_group_labels=vertical_group_labels,
                                       yticks=y_ticks, ytick_labels=y_ticks_labels.keys())
    reordered_presence_attributes = [att.replace("coverage", "present") for att in reordered_coverage_attributes]
    render_presence_heatmap_xyloglucan(merged_dataframe_alt_grouping[reordered_presence_attributes], align_cov_subtitle,
                                       name_title_insert + " Heatmap", alignment_threshold, coverage_threshold,
                                       name_prefix, graphics_output_folder,
                                       yticks=y_ticks, ytick_labels=y_ticks_labels.keys())

    render_coverage_heatmap_xyloglucan(alignment_threshold, coverage_attributes, graphics_output_folder,
                                       merged_dataframe_plot_coverage, name_prefix, name_title_insert,
                                       fixed_column_names, align_thresh_insert,
                                       f"{name_title_insert} Coverage Heatmap - gapped", alignment_threshold,
                                       vertical_group_labels=vertical_group_labels)

    # don't like this version but sort of works
    x_label = f"{name_title_insert} name"
    y_label = f"Percent with {name_title_insert} present"
    hue_label = "study_name"
    all_presence_df = merged_dataframe.reset_index()[[hue_label, any_pul_present] + reordered_presence_attributes] \
        .melt(id_vars=[hue_label], var_name=x_label, value_name=y_label)
    all_presence_df.replace(presence_column_mapper, inplace=True)
    complete_ordering = study_id_mapper.values()
    render_study_pul_proportions(merged_dataframe, count_name, graphics_output_folder, name_prefix, name_title_insert,
                                 complete_ordering, all_presence_df, "some_title",
                                 "subtitle", 0.95, 0.7, x_label, y_label, hue_label)
    # this is cleaner but broken

    render_pul_proportions(merged_dataframe, count_name, graphics_output_folder, name_prefix, name_title_insert,
                           "some_title", "subtitle", 0.95, 0.7, col_label=f"{name_title_insert} name",
                           data_label=f"Percent with {name_title_insert} present", id_label="study_name",
                           id_ordering=study_id_mapper.values(), df_colnames=reordered_presence_attributes)

    render_pul_proportions(merged_dataframe, count_name, graphics_output_folder, name_prefix, name_title_insert,
                           "some_title", "subtitle", 0.95, 0.7, col_label=f"{name_title_insert} name",
                           data_label=f"Percent with {name_title_insert} present", id_label="study_name",
                           id_ordering=study_id_mapper.values(), df_colnames=reordered_presence_attributes,
                           flip_col_id=True)

    try:
        reordered_fixed_column_names.remove('Xylanibacter\nruminicola')
        reordered_fixed_column_names.remove('Prevotella\ndentalis')
        reordered_fixed_column_names.remove('Xylanibacter\noryzae')
        reordered_fixed_column_names.remove('Prevotella\nherbatica')
        reordered_fixed_column_names.remove('Segatella\npaludivivens')
        reordered_coverage_attributes.remove('Xylanibacter ruminicola\ncoverage')
        reordered_coverage_attributes.remove('Prevotella dentalis\ncoverage')
        reordered_coverage_attributes.remove('Xylanibacter oryzae\ncoverage')
        reordered_coverage_attributes.remove('Prevotella herbatica\ncoverage')
        reordered_coverage_attributes.remove('Segatella paludivivens\ncoverage')
        reordered_presence_attributes = [att.replace("coverage", "present") for att in reordered_coverage_attributes]
    except:
        pass

    lower_cov_thresh = 0.5
    lower_cov_thresh_subtitle = f"{align_thresh_insert}, {lower_cov_thresh * 100: 2.0f}% coverage threshold"
    lower_graphics_output_folder = os.path.join(graphics_base_folder,
                                                f"coverage-{coverage_threshold}_alignment-{lower_cov_thresh}")
    if not os.path.isdir(lower_graphics_output_folder):
        os.mkdir(lower_graphics_output_folder)

    render_pul_proportions(merged_dataframe, count_name, lower_graphics_output_folder, name_prefix, name_title_insert,
                           "MLG-PUL Comparison", lower_cov_thresh_subtitle, 0.95, lower_cov_thresh,
                           col_label=f"{name_title_insert} name",
                           data_label=f"Percent with {name_title_insert} present", id_label="study_name",
                           id_ordering=study_id_mapper.values(), df_colnames=reordered_coverage_attributes,
                           legend_beside=True)

    render_pul_proportions(merged_dataframe, count_name, lower_graphics_output_folder, name_prefix, name_title_insert,
                           "MLG-PUL Comparison", lower_cov_thresh_subtitle, 0.95, lower_cov_thresh,
                           col_label=f"{name_title_insert} name",
                           data_label=f"Percent with {name_title_insert} present", id_label="study_name",
                           id_ordering=study_id_mapper.values(), df_colnames=reordered_coverage_attributes,
                           flip_col_id=True)

    merged_dataframe_extra_col = deepcopy(merged_dataframe)
    dummy_cov = "dummy_coverage"
    merged_dataframe_extra_col[dummy_cov] = 1.0
    angle = 30 if len(reordered_coverage_attributes) >= 8 else 0

    render_pul_proportions(merged_dataframe_extra_col, count_name, lower_graphics_output_folder, name_prefix, name_title_insert,
                           "MLG-PUL Comparison", lower_cov_thresh_subtitle, 0.95, lower_cov_thresh,
                           col_label=f"{name_title_insert} name",
                           data_label=f"Percent with {name_title_insert} present", id_label="study_name",
                           id_ordering=study_id_mapper.values(), df_colnames=reordered_coverage_attributes + [dummy_cov],
                           flip_col_id=True, pul_name_angle=angle, x_label="Study")

    any_colname = "Any Prevotellaceae"
    prevotella_attrs = reordered_coverage_attributes[:-1]
    bacteroides_attrs = reordered_coverage_attributes[-1]
    merged_dataframe_extra_col[any_colname] = merged_dataframe_extra_col[prevotella_attrs].max(axis=1)

    render_stacked_bar_comparison(merged_dataframe_extra_col, "MLG-PUL Comparison", lower_cov_thresh_subtitle, 0.5,
                                  "study_name", "Percent with MLG-PUL present", study_id_mapper.values(), any_colname,
                                  bacteroides_attrs, lower_graphics_output_folder)

    if coverage_comparisons:
        for col1, col2 in coverage_comparisons:
            render_coverage_comparison(merged_dataframe, col1, col2, graphics_output_folder, "All Studies")
            render_coverage_comparison_regression(merged_dataframe, col1, col2, graphics_output_folder, "All Studies")
            for study_dataframe in dataframe_dict.values():
                study_name = study_dataframe.study_name.unique()[0]
                study_output_folder = os.path.join(graphics_output_folder, study_name)
                if not os.path.exists(study_output_folder):
                    os.mkdir(study_output_folder)
                render_coverage_comparison(study_dataframe, col1, col2, study_output_folder, study_name)
                render_coverage_comparison_regression(study_dataframe, col1, col2, study_output_folder, study_name)


    # heatmaps by study
    for study_id, study_dataframe in dataframe_dict.items():
        study_name = study_dataframe.study_name.unique()[0]
        study_output_folder = os.path.join(graphics_output_folder, study_name)
        study_presence_df = study_dataframe.reset_index().set_index(heatmap_index_levels)
        study_presence_df = study_presence_df[heatmap_presence_attributes + presence_attributes]
        study_presence_df.pop(any_pul_present)
        # study_presence_df = study_presence_df.set_index("phenotype", append=True) \
        #     .reorder_levels(heatmap_index_levels) \
        #     .reindex(ordering, level=0)
        study_prefix = study_alias_dict[study_id] + ("" if name_prefix is None else name_prefix)
        if study_id == 'merged_cs_prism_metadata':
            study_presence_df = study_dataframe[stratify_properties + heatmap_presence_attributes + presence_attributes + ["antibiotics"]]
            study_presence_df = study_presence_df[study_presence_df.antibiotics == False]
            study_presence_df.pop(any_pul_present)
            study_presence_df.pop("antibiotics")
            study_presence_df = study_presence_df.set_index("phenotype", append=True).reorder_levels(heatmap_index_levels)
            study_presence_df = study_presence_df.reindex(ordering, level=0)
            study_prefix = study_alias_dict[study_id] + ", no antibiotics" + ("" if name_prefix == '' else ', ' + name_prefix)
        study_presence_df.rename(columns=presence_column_mapper, inplace=True)
        render_presence_heatmap_xyloglucan(study_presence_df, align_cov_subtitle,
                                           f"{name_title_insert} Presence Heatmap", alignment_threshold,
                                           coverage_threshold, study_prefix, study_output_folder)
        render_coverage_heatmap_xyloglucan(alignment_threshold, reordered_coverage_attributes, study_output_folder,
                                           # merged_dataframe_alt_grouping[reordered_coverage_attributes],
                                           study_dataframe[reordered_coverage_attributes],
                                           study_alias_dict[study_id], name_title_insert,
                                           reordered_fixed_column_names, align_thresh_insert,
                                           f"{name_title_insert} Coverage Heatmap", alignment_threshold,
                                           vertical_group_labels=vertical_group_labels, cell_dim=(0.03, 0.4),
                                           xtick_fontsize=4, linewidth=0.2, cbar_fraction=0.04)
        pass

    # histogram of xyGUL count
    ax = seaborn.histplot(data=merged_dataframe, x=count_name, hue="phenotype", stat="percent", common_norm=False,
                          multiple="dodge", shrink=0.8, discrete=True, hue_order=ordering, edgecolor="black")
    plt.suptitle(f"All studies count histogram")
    ax.set_title(align_cov_subtitle, fontdict={'size': 10})
    all_study_counts = phenotype_labels_with_counts(merged_dataframe, ordering)
    plt.legend(labels=list(all_study_counts.__reversed__()))
    # plt.title(f"{title_prefix}{threshold*100: 2.0f}% coverage threshold", )
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}xyGUL_count_histogram.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}xyGUL_count_histogram.png"), dpi=300)
    matplotlib.pyplot.show()

    # presence bar chart for all studies
    all_studies_presence_bargraph(all_study_counts, any_pul_present, ax, merged_dataframe, name_title_insert,
                                  ordering, presence_attributes, align_cov_subtitle, graphics_output_folder, name_prefix)

    # # treatment mean comparisons
    # treatment_group_df = merged_dataframe[["phenotype", any_pul_present] + treatments].melt(id_vars=["phenotype"] + treatments, var_name="xyGUL name", value_name="xyGUL present")
    # treatment_group_df_2 = merged_dataframe[["phenotype", any_pul_present] + treatments].dropna(how="any")
    # ax = seaborn.barplot(data=treatment_group_df_2, x=treatments, y="xyGUL present", hue="phenotype", ci=None, linewidth=1,
    #                      edgecolor='black', hue_order=ordering)
    # matplotlib.pyplot.show()

    # comparisons by study
    for study_id, study_dataframe in dataframe_dict.items():
        x_label = f"{name_title_insert} name"
        y_label = f"Percent with {name_title_insert} present"
        study_presence_df = study_dataframe[["phenotype", any_pul_present] + presence_attributes] \
            .melt(id_vars=["phenotype"], var_name=x_label, value_name=y_label)
        study_prefix = study_alias_dict[study_id] + ("" if name_prefix is None else name_prefix)
        study_presence_df.replace(presence_column_mapper, inplace=True)
        title = f"{study_prefix} {name_title_insert} comparison"
        save_study_pul_proportions(study_presence_df, graphics_output_folder, study_prefix, name_title_insert)
        render_study_pul_proportions(study_dataframe, count_name, graphics_output_folder, study_prefix,
                                     name_title_insert, ordering, study_presence_df, title, align_cov_subtitle,
                                     alignment_threshold, coverage_threshold, x_label, y_label)
        if study_id == 'merged_cs_prism_metadata':
            study_dataframe = study_dataframe[study_dataframe.antibiotics == False]
            study_presence_df = study_dataframe[["phenotype", any_pul_present] + presence_attributes] \
                .melt(id_vars=["phenotype"], var_name=x_label, value_name=y_label)
            study_prefix = study_alias_dict[study_id] + ", no antibiotics" + ("" if name_prefix == '' else ', ' + name_prefix)
            study_presence_df.replace(presence_column_mapper, inplace=True)
            title = f"{study_prefix} {name_title_insert} comparison"
            save_study_pul_proportions(study_presence_df, graphics_output_folder, study_prefix, name_title_insert)
            render_study_pul_proportions(study_dataframe, count_name, graphics_output_folder, study_prefix,
                                         name_title_insert, ordering, study_presence_df, title, align_cov_subtitle,
                                         alignment_threshold, coverage_threshold, x_label, y_label)
        elif study_id == "HMP2_metadata":
            study_dataframe = study_dataframe.query("GigaBases > 1 and antibiotics == False")
            study_presence_df = study_dataframe[["phenotype", any_pul_present] + presence_attributes] \
                .melt(id_vars=["phenotype"], var_name=x_label, value_name=y_label)
            study_prefix = study_alias_dict[study_id] + ", no antibiotics, 1 gigabase minimum" + ("" if name_prefix == '' else ', ' + name_prefix)
            study_presence_df.replace(presence_column_mapper, inplace=True)
            title = f"{study_prefix} {name_title_insert} comparison"
            save_study_pul_proportions(study_presence_df, graphics_output_folder, study_prefix, name_title_insert)
            render_study_pul_proportions(study_dataframe, count_name, graphics_output_folder, study_prefix,
                                         name_title_insert, ordering, study_presence_df, title, align_cov_subtitle,
                                         alignment_threshold, coverage_threshold, x_label, y_label)

        # treatment comparison charts
        for treatment in treatments:
            try:
                treatment_presence_df = merged_dataframe[["phenotype", treatment, any_pul_present] + presence_attributes]
                treatment_presence_df = treatment_presence_df[pandas.notnull(treatment_presence_df[treatment])]
                if len(treatment_presence_df[(treatment_presence_df[treatment] == True)]) > 0:
                    plot_treatment_comparison(ax, graphics_output_folder, name_prefix, name_title_insert, ordering,
                                              align_cov_subtitle, treatment, treatment_presence_df)
            except KeyError:
                print(f"{treatment} not in {study_id} data")
            except IndexError:
                print(f"{treatment} IndexError in {study_id} data. Probably no patients who got this treatment.")
    print("Done rendering!")  # todo add back breakpoint


def phenotype_labels_with_counts(dataframe, ordering):
    counts = dataframe.groupby("phenotype").count().reindex(ordering)["Bases"]
    count_label = [f"{patient_category:<8} n ={n:>3}" for patient_category, n in counts.items()]
    return count_label


def render_bases_coverage_by_study(dataframe, cazyme_name, graphics_output_folder, name_prefix, study_dict, title,
                                   subtitle):
    # # "Melt" the dataset to "long-form" or "tidy" representation
    # iris = pd.melt(iris, "species", var_name="measurement")
    dataframe = dataframe.replace(study_dict)

    # Initialize the figure
    fig, ax = plt.subplots()
    # seaborn.despine(bottom=True, left=True)

    # Show each observation with a scatterplot
    seaborn.stripplot(
        data=dataframe, x="GigaBases", y="study_name", hue="phenotype",
        dodge=True, alpha=.25, zorder=1, ax=ax, legend=False
    )
    # handles, labels = ax.get_legend_handles_labels()

    # Show the conditional means, aligning each pointplot in the
    # center of the strips by adjusting the width allotted to each
    # category (.8 by default) by the number of hue levels
    seaborn.pointplot(
        data=dataframe, x="GigaBases", y="study_name", hue="phenotype",
        join=False, dodge=.8 - .8 / 3,
        markers="d", errorbar=None, ax=ax
    )

    ax.figure.suptitle(f"{name_prefix} Metagenome Sequencing Depth by Study and Phenotype")
    matplotlib.pyplot.title(subtitle)
    ax.set_ylabel("Study")
    ax.set_xlabel("Base pairs (GigaBases)")
    ax.legend().set_title("Phenotype")

    # # Improve the legend
    # ax.get_legend().remove()
    fig.tight_layout()
    fig.set_dpi(200)
    # fig.legend(handles, labels, title="Phenotype", loc='upper right', ncol=3, bbox_to_anchor=(.75, 0.98))
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}sequencing_depth_all_studies.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}sequencing_depth_all_studies.png"), dpi=300)
    matplotlib.pyplot.show()
    pass


def render_bases_coverage_scatterplot(dataframe, cazyme, graphics_output_folder, name_prefix, title, subtitle,
                                      ordering):
    cazyme_insert = 'xygul' if cazyme == 'genome' or cazyme is None or cazyme == '' else cazyme
    y_var = f"max_{cazyme_insert}_coverage"
    grid = seaborn.lmplot(data=dataframe, x="GigaBases", y=y_var, hue='phenotype', hue_order=ordering, logistic=True,
                          aspect=1.3, legend=True)
    grid.ax.set_xlim(left=-0.1)
    grid.figure.set_dpi(200)
    grid.figure.suptitle(f"{name_prefix} {title}")
    matplotlib.pyplot.title(subtitle)
    grid.facet_axis(0, 0).set_ylabel(f"Highest coverage from any {cazyme_insert}")
    grid.facet_axis(0, 0).set_xlabel("Metagenome sample base pair count (GigaBases)")
    grid.legend.set_frame_on(True)
    grid.legend.set_title("Phenotype")
    count_labels = phenotype_labels_with_counts(dataframe, ordering)
    for i, text_obj in enumerate(grid.legend.texts):
        text_obj._text = count_labels[i]
    seaborn.despine(right=False, top=False)
    grid.ax.axvline(1, color="black", linestyle="--", linewidth=1.2)
    grid.figure.tight_layout(rect=[0.03, 0.03, 1, 1])  # left x, bottom y, right x, top y
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}base_coverage_scatter_all_studies.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}base_coverage_scatter_all_studies.png"), dpi=300)
    matplotlib.pyplot.show()


def plot_treatment_comparison(ax, graphics_output_folder, name_prefix, name_title_insert, ordering, subtitle, treatment,
                              treatment_presence_df):
    seaborn.set(font="monospace")

    true_counts = treatment_presence_df.query(f"{treatment} == True").groupby("phenotype").count().reindex(ordering)[
        treatment]
    true_counts = true_counts.fillna(0).astype(int)
    false_counts = treatment_presence_df.query(f"{treatment} == False").groupby("phenotype").count().reindex(ordering)[
        treatment]
    false_counts = false_counts.fillna(0).astype(int)
    true_counts_labels = [f"{patient_category:<8} n ={n:>3}" for patient_category, n in true_counts.items()]
    false_counts_labels = [f"{patient_category:<8} n ={n:>3}" for patient_category, n in false_counts.items()]
    treatment_presence_df = treatment_presence_df.melt(id_vars=["phenotype", treatment],
                                                       var_name=f"{name_title_insert} name",
                                                       value_name=f"Percent with {name_title_insert} present")
    # todo : fix ylimits to be [0,100.5] eg. ax.set_ylim([0, 100.5])
    grid = seaborn.FacetGrid(treatment_presence_df, col=treatment)
    grid.map_dataframe(seaborn.barplot, x=f"{name_title_insert} name", y=f"Percent with {name_title_insert} present",
                       hue="phenotype", errorbar=None, linewidth=1,
                       edgecolor='black', palette="deep", hue_order=ordering,
                       estimator=lambda x: np.mean(x) * 100)
    # grid.add_legend()
    # for ax in grid.axes.ravel():
    grid.facet_axis(0, 0).legend(labels=false_counts_labels, prop={'size': 5})
    grid.facet_axis(0, 1).legend(labels=true_counts_labels, loc="upper right", prop={'size': 5})

    grid.facet_axis(0, 0).set_title(grid.facet_axis(0, 0).get_title(), {"fontsize": "small"})
    grid.facet_axis(0, 1).set_title(grid.facet_axis(0, 1).get_title(), {"fontsize": "small"})
    # ax = seaborn.barplot(data=treatment_presence_df, x="xyGUL name", y="xyGUL present", hue="phenotype", ci=None, linewidth=1, edgecolor='black')
    grid.set_xticklabels(grid.facet_axis(0, 0).get_xticklabels(), rotation=0, ha="center", fontsize=3)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=5, ha="right", fontsize=7)
    grid.figure.set_dpi(300)
    grid.figure.tight_layout()

    # grid.figure.suptitle(f"{name_title_insert} {treatment} comparison", y=0.99)
    grid.figure.suptitle(subtitle, y=0.99, size="small")
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}treatment_comparison_{treatment}.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}treatment_comparison_{treatment}.png"), dpi=300)
    matplotlib.pyplot.show()


def all_studies_presence_bargraph(all_study_counts, any_pul_present, ax, merged_dataframe, name_title_insert, ordering,
                                  presence_attributes, subtitle, graphics_output_folder, name_prefix,
                                  x_name="xyGUL name", y_name="xyGUL present", hue_name="phenotype"):
    presence_df = merged_dataframe[["phenotype", any_pul_present] + presence_attributes].melt(id_vars=[hue_name],
                                                                                              var_name=x_name,
                                                                                              value_name=y_name)
    ax = seaborn.barplot(data=presence_df, x=x_name, y=y_name, hue=hue_name, errorbar=None, linewidth=1,
                         edgecolor='black', hue_order=ordering, estimator=lambda x: np.mean(x) * 100)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha="center", fontsize=7)
    ax.set_ylim([0, 100.5])
    plt.legend(labels=all_study_counts)
    ax.figure.suptitle(f"All studies {name_title_insert} presence")
    pul_name_angle, x_labels = adjusted_labels(ax.get_xticklabels(), name_title_insert)
    ax.set_xticklabels(x_labels, rotation=pul_name_angle, fontsize=7)
    ax.set_ylabel("Percent with XyGUL present")
    ax.set_title(subtitle, fontdict={'size': 10})
    ax.figure.set_dpi(300)
    ax.figure.tight_layout()
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}{name_title_insert}_presence_all_studies.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{name_prefix}{name_title_insert}_presence_all_studies.png"), dpi=300)
    matplotlib.pyplot.show()


def save_study_pul_proportions(study_presence_df, graphics_output_folder, name_prefix, cazyme):
    result = study_presence_df.groupby([f"{cazyme} name", "phenotype"]).mean()
    print(name_prefix, f"{cazyme} comparison")
    print(result)
    result.to_csv(os.path.join(graphics_output_folder, f"{name_prefix}_xygul_comparison.tsv"), sep='\t')


def render_study_pul_proportions(study_dataframe, count_name, graphics_output_folder, name_prefix, name_title_insert,
                                 hue_ordering, study_presence_df, title, subtitle, perc_ident_threshold,
                                 coverage_threshold, xlabel, ylabel, hue_label="phenotype"):

    ax = seaborn.barplot(data=study_presence_df, x=xlabel, y=ylabel, hue=hue_label,
                         errorbar=None, hue_order=hue_ordering, edgecolor="#000000",
                         estimator=lambda x: np.mean(x) * 100)

    study_counts = study_dataframe.groupby(hue_label).count().reindex(hue_ordering)[count_name]
    study_counts = study_counts.fillna(0).astype(int)
    study_counts_labels = [f"{patient_category:<8} n ={n:>3}" for patient_category, n in study_counts.items()]
    # oldtitle = f"{f'{name_prefix} ' if name_prefix else ''}{name_title_insert} comparison"
    ax.figure.suptitle(title)
    plt.legend(labels=study_counts_labels, loc="upper right")
    plt.title(subtitle)
    pul_name_angle, x_labels = adjusted_labels(ax.get_xticklabels(), name_title_insert)
    ax.set_xticklabels(x_labels, rotation=pul_name_angle, ha="center", fontsize=7)
    ax.set_ylim([0, 100.5])
    ax.figure.set_dpi(300)
    ax.figure.tight_layout()
    ax.set_position([0.1266, 0.1522, 0.971875-0.1266, 0.839583-0.1522])
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.png"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.eps"), dpi=300)
    plt.show()
    pass


def render_pul_proportions(dataframe, count_name, graphics_output_folder, name_prefix, name_title_insert, title,
                           subtitle, perc_ident_threshold, coverage_threshold, col_label, data_label, id_label,
                           id_ordering, df_colnames, flip_col_id=False, legend_beside=False, font_size=7,
                           pul_name_angle=0, x_label=None):

    all_presence_df = dataframe.reset_index()[[id_label] + list(df_colnames)] \
        .melt(id_vars=[id_label], var_name=col_label, value_name=data_label)
    # apply threshold to coverage data, doesn't effect already thresholded presence data
    all_presence_df[data_label] = all_presence_df[data_label].apply(lambda num: 1 if num >= coverage_threshold else 0)
    # all_presence_df.replace(presence_column_mapper, inplace=True)

    # complete_ordering = study_id_mapper.values()
    if flip_col_id:
        col_label, id_label = id_label, col_label
        id_ordering, df_colnames = df_colnames, id_ordering

    if x_label:
        all_presence_df.rename(columns={col_label: x_label}, inplace=True)
        col_label = x_label

    ax = seaborn.barplot(data=all_presence_df, x=col_label, y=data_label, hue=id_label,
                         errorbar=None, hue_order=id_ordering, edgecolor="#000000",
                         estimator=lambda x: np.mean(x) * 100)

    if not flip_col_id:
        study_counts = dataframe.groupby(id_label).count().reindex(id_ordering)[count_name]
        study_counts = study_counts.fillna(0).astype(int)
        legend_labels = [f"{patient_category:<8} n ={n:>3}" for patient_category, n in study_counts.items()]
    else:
        legend_labels = [string.replace("\npresent", '') for string in id_ordering]
        legend_labels = [string.replace("\ncoverage", '') for string in legend_labels]
    # oldtitle = f"{f'{name_prefix} ' if name_prefix else ''}{name_title_insert} comparison"
    ax.figure.suptitle(title)
    if legend_beside:
        plt.legend(labels=legend_labels, loc="upper right", bbox_to_anchor=(1.265, 1), fontsize=font_size)
    else:
        plt.legend(labels=legend_labels, loc="upper right", fontsize=font_size)
    plt.title(subtitle)
    # pul_name_angle, x_labels = adjusted_labels(ax.get_xticklabels(), name_title_insert)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=pul_name_angle, ha="center", fontsize=font_size)
    ax.set_ylim([0, 100.5])
    ax.figure.set_dpi(300)
    ax.figure.tight_layout()

    if legend_beside:
        ax.set_position([0.1266, 0.1522, 0.971875-0.28, 0.839583-0.1522])
    else:
        ax.set_position([0.1266, 0.1522, 0.971875-0.1266, 0.839583-0.1522])
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.png"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.pdf"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{title} {perc_ident_threshold}_perc_ident {coverage_threshold}_cov.svg"), dpi=300)
    plt.show()
    pass


def render_xyloglucan_timeseries_graphs(dataframes, coverage_threshold, pul_names, alignment_threshold, cazyme,
                                        study_names_dict, graphics_output_folder, name_prefix="",
                                        presence_attributes=None, coverage_attributes=None):
    ordering = ["Control", "UC", "CD", "NA"]
    if not os.path.exists(graphics_output_folder):
        os.mkdir(graphics_output_folder)
    cazyme_insert = 'xygul' if cazyme == 'genome' or cazyme is None or cazyme == '' else cazyme
    name_title_insert = "XyGUL" if cazyme == 'genome' or cazyme is None or cazyme == '' else cazyme
    max_variable = f"max_{cazyme_insert}_coverage"
    any_xygul_attribute = f"Any\n{'XyGUL' if cazyme == 'genome' else cazyme}\npresent"
    pul_names.append(max_variable)
    align_thresh_insert = f"{alignment_threshold * 100: 2.0f}% percent identity threshold"
    coverage_thresh_insert = f"{coverage_threshold * 100: 2.0f}% coverage threshold"
    align_cov_subtitle = f"{align_thresh_insert}, {coverage_thresh_insert}"

    for study_name, study_data in dataframes.items():
        x_label = f"{name_title_insert} name"
        y_label = f"Percent with {name_title_insert} present"
        study_presence_by_sample_df = study_data[["phenotype", any_xygul_attribute] + presence_attributes]\
            .melt(id_vars=["phenotype"], var_name=x_label, value_name=y_label)
        patient_codes = {}
        title = f"{study_name} by sample {name_title_insert} comparison"
        study_prefix = f"{study_name} by sample"
        #todo: remove line: save_study_pul_proportions(study_data, graphics_output_folder, study_prefix, name_title_insert)
        save_study_pul_proportions(study_presence_by_sample_df, graphics_output_folder, study_prefix, name_title_insert)
        render_study_pul_proportions(study_data, "study_name", graphics_output_folder, study_prefix,
                                     name_title_insert, ordering, study_presence_by_sample_df, title,
                                     align_cov_subtitle, alignment_threshold, coverage_threshold, x_label, y_label)
        if study_name == "HMP2_metadata":
            hmp2_patient_notes = {"C3020": "Only one sample, which had low sequencing depth (~0.85 Gigabases).",
                                  "C3033": "Only one sample, which had low sequencing depth (~0.1 Gigabases).",
                                  "H4039": "Child of unknown age (age < 13 years). 12 visits, good sequencing depth.",
                                  "M2010": "Only three samples, has max GH5_4 coverage of ~33% despite being below full "
                                           "XyGUL coverage",
                                  "P6024": "Adolescent ( 13 <= age < 18 years), about half of samples have low "
                                           "sequencing depth(< 1 Gigabases)",
                                  "M2034": "Adult female, 20 visits, reasonable sequencing depth, no antibiotics, "
                                           "no obvious explanation for lack of XyGUL over most of study.",
                                  "C3002": "Senior (65+) CD patient, sequencing depth looks fine, no antibiotics "
                                           "or immunosuppressants, 15 visits. No obvious explanation for lack of XyGUL "
                                           "over most of study.",
                                  "C3008": "CD patient, received antibiotics",
                                  "C3009": "Adult male, 12 visits, no antibiotics, sequencing depth looks fine. "
                                           "No obvious explanation for lack of XyGUL over most of study.",
                                  "C3030": "Adult male CD patient, 10 visits, no antibiotics, no obvious explanation.",
                                  "M2083": "Adult Male UC patient. 17 visits. No explanation.",
                                  "M2069": "UC patient, received antibiotics",
                                  "P6018": "Control adolescent ( 13 <= age < 18 years), 24 visits, good sequencing depth. "
                                           "No explanation.",
                                  "M2072": "Control Adult Male, 24 visits, good sequencing depth. No explanation.",
                                  "H4038": "CD Adolescent ( 13 <= age < 18 years), has ~35% XyGUL coverage, but near "
                                           "zero GH5_4 coverage.",
                                  "H4032": "CD child, has ~35% XyGUL coverage, but near "
                                           "zero GH5_4 coverage."
                                  }
            # patient codes: s = low sample count
            #                l = low sequencing depth
            #                c = child
            #                a = antibiotic treatment
            #                x = XyGUL present but no GH5_4
            hmp2_patient_codes = {
                "C3020": "s,l",
                "C3033": "s,l",
                "H4039": "c",
                "M2010": "s",
                "C3008": "a",
                "M2069": "a",
                "H4038": "x",
                "H4032": "x"
            }
            patient_codes = hmp2_patient_codes
            study_data["Visit"] = map_df_timeseries_count(study_data,
                                                          patient_id_varname="patient_id",
                                                          timepoint_varname="timepoint")
            study_data["Approximate biweekly sample"] = rescale_sparse_timeseries_indices(study_data, "patient_id", "timepoint")
            timepoint_name = "Approximate biweekly sample"
            xlabel = "Biweekly sample"
        else:
            timepoint_name = "timepoint"
            xlabel = None

        for attribute in presence_attributes:
            heatmap_data_old = study_data[[timepoint_name, "phenotype", "patient_id", attribute]]\
                            .pivot(columns=timepoint_name,
                                   index=["phenotype", "patient_id"],
                                   values=attribute)
            heatmap_data_old["Any"] = heatmap_data_old.max(axis=1)
            heatmap_data_old = heatmap_data_old.reindex(ordering, level=0)
            heatmap_data = format_heatmap_data(study_data, coverage_threshold, timepoint_name, attribute, ordering,
                                               max_variable)
            # todo: use the timepoint name from the tsv file, whether it's month, week, timepoint, visit number, etc
            species_name = ' '.join(attribute.split('\n')[0:-1])
            if species_name == "Any XyGUL":
                species_name = "Any"
            render_presence_heatmap_xyloglucan(heatmap_data, align_cov_subtitle,
                                               f"{study_names_dict[study_name]} {cazyme_insert} timeseries", alignment_threshold,
                                               coverage_threshold, name_prefix, graphics_output_folder, is_binary=True,
                                               xlabel=xlabel, species_name=species_name, no_data_colour="white",
                                               row_annotations=patient_codes, cbar_label=f"{species_name} {cazyme_insert} presence")

        patients_ever_below_thresh = pd.Series(study_data[["patient_id"] + coverage_attributes]
                                               .groupby("patient_id")
                                               .min()
                                               .query(f"max_{cazyme}_coverage < {coverage_threshold}")
                                               .index)
        study_data_subset = study_data.merge(patients_ever_below_thresh, on="patient_id")
        for name in pul_names:
            render_timeseries_lineplot(cazyme_insert, graphics_output_folder, max_variable, name, name_prefix,
                                       study_data)


def format_heatmap_data(study_data, coverage_threshold, timepoint_name, attribute, ordering, max_variable):
    if attribute.__contains__("Any"):
        cov_attribute = max_variable
    else:
        cov_attribute = attribute.replace("present", "coverage")
    heatmap_data = study_data[[timepoint_name, "phenotype", "patient_id", cov_attribute]] \
        .pivot(columns=timepoint_name,
               index=["phenotype", "patient_id"],
               values=cov_attribute)
    heatmap_data["Any"] = heatmap_data.max(axis=1)
    heatmap_data = heatmap_data.reindex(ordering, level=0)
    return heatmap_data.applymap(lambda x: numpy.nan if pandas.isna(x) else 1 if x > coverage_threshold else 0)


def rescale_sparse_timeseries_indices(timeseries_dataframe: pandas.DataFrame, patient_id_varname: str, timepoint_varname: str,
                                      min_index: int = None, max_index: int = None):
    sample_counts = map_df_timeseries_count(timeseries_dataframe,
                                            patient_id_varname=patient_id_varname,
                                            timepoint_varname=timepoint_varname)
    scaling_ratio = max(sample_counts) / max(timeseries_dataframe[timepoint_varname])
    rescaled_timepoints = round(timeseries_dataframe[timepoint_varname] * scaling_ratio).astype(int)
    sample_min = min(sample_counts)
    sample_max = max(sample_counts)
    if not min_index:
        min_index = sample_min
    if not max_index:
        max_index = sample_max
    if max_index - min_index > sample_max - sample_min:
        raise Exception("Desired reindexing is narrower than minimum data width of your timeseries!\n"
                        "At least one patient (i.e. row) has more entries than will fit in the desired reindexing.")
    return diffuse_indices(rescaled_timepoints, min_index, max_index)


def map_df_timeseries_count(study_data, patient_id_varname, timepoint_varname):
    return numpy.concatenate(
        list(map(lambda x: x - min(x) + 1, study_data.groupby(patient_id_varname)[timepoint_varname].indices.values()))
    )


def diffuse_indices(array: pandas.Series, min_index=None, max_index=None, inplace=False):
    if inplace:
        array_copy = array
    else:
        array_copy = array.copy()

    if not min_index:
        min_index = min(array_copy)
    else:
        array_copy[array_copy < min_index] = min_index

    if not max_index:
        max_index = max(array_copy)
    else:
        array_copy[array_copy > max_index] = max_index

    has_duplicates = True

    iteration = 0
    while has_duplicates:
        has_duplicates = False
        if iteration % 2 == 0:
            for i in range(len(array_copy)):
                if i < len(array_copy)-1 and array_copy[i] == array_copy[i+1]:
                    has_duplicates = True
                    if array_copy[i+1] < max_index:
                        j = i+1
                        # need to loop ahead using j to increment multiples indices when 3 or more indices
                        # are equal at once
                        while j < len(array_copy) and array_copy[j] == array_copy[i]:
                            array_copy[j] += 1
                            j += 1
        else:
            for i in range(len(array_copy)-1, -1, -1):
                if i > 0 and array_copy[i] == array_copy[i-1]:
                    has_duplicates = True
                    if array_copy[i-1] > min_index:
                        j = i-1
                        while j >= 0 and array_copy[j] == array_copy[i]:
                            array_copy[j] -= 1
                            j -= 1
        iteration += 1

    return array_copy


def render_timeseries_lineplot(cazyme_insert, graphics_output_folder, max_variable, xygul_name, name_prefix,
                               study_data):
    y_label = f"{xygul_name}\n{cazyme_insert}\ncoverage" if xygul_name != max_variable else max_variable
    spaced_xygul_name = xygul_name.replace('\n', ' ')
    study_name = study_data["study_name"].unique()[0]
    ax = seaborn.lineplot(data=study_data, x="timepoint", y=y_label, hue="patient_id")
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.figure.tight_layout(rect=[0, 0.03, 1, 0.95])
    ax.set_xticks(range(1, 13))
    ax.set_ylim([0, 1])
    ax.set_xlabel("Month")
    if y_label == max_variable:
        ax.set_ylabel("Highest coverage of any XyGUL")
    ax.set_title(f"{study_name} timeseries, 1 month intervals")
    ax.figure.set_dpi(300)
    filename = f"{name_prefix} {study_name} {spaced_xygul_name}_{cazyme_insert}_coverage_timeseries"
    plt.savefig(os.path.join(graphics_output_folder, f"{filename}.png"), dpi=300)
    plt.savefig(os.path.join(graphics_output_folder, f"{filename}.pdf"), dpi=300)
    matplotlib.pyplot.show()


def analyze_multiple_perc_ident_levels(levels, coverage_threshold, metadata_files, pul_dataframe_file, pul_list, pul_id,
                                       cazymes_to_analyze, pul_alias_dict, study_name_alias_dict, study_id_dict,
                                       pul_names, mixed_loci_name, computed_df_folder, graphics_out_folder):
    dataframes = {}
    aggregated_df = pandas.DataFrame()
    files_failed = 0
    for level in levels:
        data_folder = os.path.join(os.getcwd(), "xyloglucan output", f"{level}")
        dataframes[level], temp_files_failed, pul_list, coverage_attributes, presence_attributes = \
            compute_xyloglucan_data_df(metadata_files, computed_df_folder, pul_list, pul_dataframe_file, pul_id,
                                       data_folder=data_folder, coverage_threshold=coverage_threshold,
                                       cazyme_list=cazymes_to_analyze, pul_alias_dict=pul_alias_dict,
                                       mixed_loci_name=mixed_loci_name)
        files_failed += temp_files_failed
        # if len(dataframes[level]) == 1:
        for key in dataframes[level].keys():
            # metadata = metadata[key]
            aggregated_df[f"Any\n{pul_id}\nPresent\n{level} Perc. Ident."] = dataframes[level][key][f"Any\n{pul_id}\npresent"]
            aggregated_df[f"Any\n{pul_id}\nPresent\n{level} Perc. Ident."] = dataframes[level][key][f"Any\n{pul_id}\npresent"]

        for cazyme in cazymes_to_analyze:
            graph_rendering_pipeline(dataframes[level], study_name_alias_dict, study_id_dict, pul_list, pul_names,
                                     coverage_attributes[cazyme], presence_attributes[cazyme], cazyme,
                                     coverage_threshold, float(level) / 100, graphics_out_folder, cazyme, False, False)

    print("Done multiple percent identity level comparison!")


def analyze_xyloglucan_heatmap():
    metadata_files = ["xyloglucan_metadata\\PRJEB2054_SraRunTable.txt",
                      "xyloglucan_metadata\\PRJNA32089_WGS_SraRunTable.csv",
                      "xyloglucan_metadata\\PRJNA175224_SraRunTable.csv",
                      "xyloglucan_metadata\\PRJNA400072_SraRunTable.csv"
                      ]
    study_names = ["MetaHIT", "Obese-Lean twins", "2012 IBD dysfunction", "CS-PRISM from HMBR"
                   ]
    pul_list = ["B_cellulosilyticus_DSM_14838_xyGUL", "B_ovatus_ATCC_8483_xyGUL",
                "B_uniformis_ATCC_8492_xyGUL_1", "B_uniformis_ATCC_8492_xyGUL_2"]
    pul_names = ["B. cellulosilyticus\nDSM 14838\nxyGUL", "B. ovatus \nATCC 8483\nxyGUL",
                 "B. uniformis\nATCC 8492\nxyGUL 1", "B. uniformis\nATCC 8492\nxyGUL 2"]
    data_folder = os.path.join(os.getcwd(), "xyloglucan output")
    threshold = 0.7
    xyloglucan_presence_data, study_metadata, files_failed = compute_xyloglucan_data(metadata_files, pul_list, data_folder=data_folder, threshold=threshold)
    print(f"WARNING: {files_failed} data files could not be read!")
    print("Analysis Done, graphing")
    title = f"{threshold*100: 2.0f}% threshold XyGUL Heatmap"
    render_heatmap_xyloglucan(xyloglucan_presence_data, study_metadata, study_names, pul_list, pul_names, title)


def analyze_hmp2_xyloglucan_puls():
    study_metadata_files = [str(project_folder / "Data" / "Metagenomics" / "HMP2Metadata" / "HMP2_metadata.tsv")]
    timeseries_metadata_files = study_metadata_files
    study_name_alias_dict = ibd_study_name_alias_dict

    study_id_dict = ibd_study_id_dict

    pul_list = ["B_ovatus_ATCC_8483_xyGUL",
                "B_cellulosilyticus_DSM_14838_xyGUL",
                "B_uniformis_ATCC_8492_xyGUL_1",
                "B_uniformis_ATCC_8492_xyGUL_2",
                "B_fluxus_YIT12057_xyGUL",
                ]
    pul_names = ["B. ovatus \nATCC 8483\nxyGUL",
                 "B. cellulosilyticus\nDSM 14838\nxyGUL",
                 "B. uniformis\nATCC 8492\nxyGUL 1",
                 "B. uniformis\nATCC 8492\nxyGUL 2",
                 "B. fluxus\nYIT12057\nXyGUL",
                 ]
    pul_alias_dict = {"B_cellulosilyticus_DSM_14838_xyGUL": "B.\ncellulosilyticus",
                      "B_ovatus_ATCC_8483_xyGUL": "B. ovatus",
                      "B_uniformis_ATCC_8492_xyGUL_1": "B. uniformis\nxyGUL 1",
                      "B_uniformis_ATCC_8492_xyGUL_2": "B. uniformis\nxyGUL 2",
                      "B_fluxus_YIT12057_xyGUL": "B. fluxus"
                      }
    pul_dataframe_file = str(project_folder / "Data" / "Metagenomics" / "XyGUL metadata (glycocage paper set).tsv")
    pickle_folder = str(project_folder / "Pickles")
    cazymes_to_analyze = ("GH5_4",)
    mixed_loci_name = "8_xyguls"
    pul_tag = "XyGUL"

    data_folder_95 = str(project_folder / "Data" / "Metagenomics" / "SampleXyGULCoverage")

    graphics_output_folder = project_folder / "Code" / "GraphOutput"
    graphics_folder_70 = str(graphics_output_folder / "70PercentCoverage")
    graphics_folder_20 = str(graphics_output_folder / "20PercentCoverage")
    graphics_folder_timeseries = str(graphics_output_folder / "timeseries")

    strict_coverage_threshold = 0.7
    loose_coverage_threshold = 0.2
    timeseries_dataframes, timeseries_files_failed, timeseries_pul_list, timeseries_coverage_attributes, timeseries_presence_attributes = \
        compute_xyloglucan_data_df(timeseries_metadata_files, pickle_folder, pul_list, pul_dataframe_file, pul_tag,
                                   data_folder=data_folder_95, coverage_threshold=loose_coverage_threshold,
                                   cazyme_list=cazymes_to_analyze, pul_alias_dict=pul_alias_dict,
                                   keep_patient_timeseries=True, mixed_loci_name=mixed_loci_name)

    # analyze_multiple_perc_ident_levels(perc_ident_levels, 0.7, study_metadata_files, pul_dataframe_file, pul_list,
    #                                    pul_tag, cazymes_to_analyze, pul_alias_dict, study_name_alias_dict,
    #                                    study_id_dict, pul_names, mixed_loci_name, pickle_folder, graphics_output_folder)
    #
    # analyze_multiple_perc_ident_levels(perc_ident_levels, 0.2, study_metadata_files, pul_dataframe_file, pul_list,
    #                                    pul_tag, cazymes_to_analyze, pul_alias_dict, study_name_alias_dict,
    #                                    study_id_dict, pul_names, mixed_loci_name, pickle_folder, graphics_output_folder)

    study_dataframes, files_failed, pul_list, coverage_attributes, presence_attributes = \
        compute_xyloglucan_data_df(study_metadata_files, pickle_folder, pul_list, pul_dataframe_file, "XyGUL",
                                   data_folder=data_folder_95, coverage_threshold=strict_coverage_threshold,
                                   cazyme_list=cazymes_to_analyze, pul_alias_dict=pul_alias_dict,
                                   mixed_loci_name=mixed_loci_name)

    print(f"WARNING: {files_failed} data files could not be read!")
    # print("Analysis Done, graphing")

    render_xyloglucan_timeseries_graphs(timeseries_dataframes, loose_coverage_threshold,
                                        [pul_alias_dict[item] for item in pul_list], alignment_threshold=0.95,
                                        cazyme="GH5_4", study_names_dict=study_name_alias_dict,
                                        graphics_output_folder=graphics_folder_timeseries,
                                        presence_attributes=timeseries_presence_attributes["GH5_4"],
                                        coverage_attributes=timeseries_coverage_attributes["GH5_4"])

    graph_rendering_pipeline(study_dataframes, study_name_alias_dict, study_id_dict, pul_list, pul_names,
                             coverage_attributes["GH5_4"], presence_attributes["GH5_4"], "GH5_4", strict_coverage_threshold,
                             0.95, graphics_folder_70, pul_tag, False, False,
                             study_id_mapper=study_id_dict)

    graph_rendering_pipeline(study_dataframes, study_name_alias_dict, study_id_dict, pul_list, pul_names,
                             coverage_attributes["GH5_4"], presence_attributes["GH5_4"], "GH5_4", loose_coverage_threshold,
                             0.95, graphics_folder_20, pul_tag, False, False,
                             study_id_mapper=study_id_dict)


    pass

if __name__ == "__main__":
    analyze_hmp2_xyloglucan_puls()
