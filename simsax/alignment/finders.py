from Bio import pairwise2

from simsax.alignment.model import AlignmentMotif, Alignment
from simsax.alignment.score_functions import create_score_functions

import numpy as np


def create_sax_signature(sax):
    return "%d_%d_%d" % (sax.window, sax.nbins, sax.nlevels)


def calculate_global_project_alignment(project_sequence_a, project_sequence_b,
                                       sax_signature, min_window_dist=0, similarity=0.95):
    # check if both sequences have calculated windows for a signature
    if sax_signature not in project_sequence_a.project_dna_sax.keys() or \
                    sax_signature not in project_sequence_a.project_dna_sax.keys():
        raise Exception("Calculate windows for a given sax signature")

    sax = project_sequence_a.project_dna_sax[sax_signature]

    score_func, gap_func = create_score_functions(sax)

    windows_a = project_sequence_a.project_dna[sax_signature]
    windows_b = project_sequence_b.project_dna[sax_signature]

    result = list()

    threshold = sax.nbins * similarity

    comparisons = len(windows_a) * len(windows_b)
    inform_level = comparisons // 25 if comparisons // 25 > 0 else 1
    progress = 0

    for i, window_a in enumerate(windows_a, 1):
        for j, window_b in enumerate(windows_b):
            if abs(window_b.start_index - window_a.start_index) >= min_window_dist or \
                            project_sequence_a.project != project_sequence_b.project:

                # print(i, j)
                alignments = pairwise2.align.globalcc(window_a.word, window_b.word,
                                                      score_func, gap_func, gap_func,
                                                      score_only=False, one_alignment_only=True)
                max_score = alignments[0][2]
                similarity_perc_score = max_score / sax.nbins
                if max_score >= threshold:
                    alignment = Alignment(window_a, window_b, alignments[0], similarity_perc_score)
                    result.append(alignment)
            progress += 1
            if progress % inform_level == 0 or progress == comparisons:
                print("{}/{}".format(progress, comparisons))
    return result


def calculate_project_alignment(project_sequence_a, project_sequence_b,
                                       sax_signature, min_window_dist=0, similarity=1, verbose=False):
    # check if both sequences have calculated windows for a signature
    if sax_signature not in project_sequence_a.project_dna_sax.keys() or \
                    sax_signature not in project_sequence_a.project_dna_sax.keys():
        raise Exception("Calculate windows for a given sax signature")

    sax = project_sequence_a.project_dna_sax[sax_signature]

    windows_a = project_sequence_a.project_dna[sax_signature]
    windows_b = project_sequence_b.project_dna[sax_signature]

    result = list()

    threshold = sax.nbins * similarity

    comparisons = len(windows_a) * len(windows_b)
    inform_level = comparisons // 25 if comparisons // 25 > 0 else 1
    progress = 0

    for i, window_a in enumerate(windows_a, 1):
        for j, window_b in enumerate(windows_b):

            if abs(window_b.start_index - window_a.start_index) >= min_window_dist or \
                            project_sequence_a.project != project_sequence_b.project:

                max_score = np.sum([1 if a == b else 0 for a, b in list(zip(list(window_a.word), list(window_b.word)))])
                similarity_perc_score = max_score / sax.nbins
            
                if max_score >= threshold:
                    alignment_tuple = (window_a.word, window_b.word, max_score, 0, len(window_a.word))
                    alignment = Alignment(window_a, window_b, alignment_tuple, similarity_perc_score)
                    result.append(alignment)
            progress += 1
            if verbose and (progress % inform_level == 0 or progress == comparisons):
                print("{}/{}".format(progress, comparisons))
    return result


def find_motifs(input_alignments):
    motifs = []

    for align in input_alignments:
        found_motif = False
        for motif in motifs:
            if (motif.alignment_belongs_to_motif(align) == False) and \
                    (motif.word_belongs_to_motif(align.window_a.word) or motif.word_belongs_to_motif(
                        align.window_b.word)):
                motif.add_alignment(align)
                found_motif = True

        if found_motif == False:
            new_motif = AlignmentMotif()
            new_motif.add_alignment(align)
            motifs.append(new_motif)

    return motifs

