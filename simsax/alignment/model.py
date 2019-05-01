from copy import copy
import pandas as pd


class ProjectSequence(object):
    def __init__(self, project, get_sequence):
        self.project = project
        self.sequence = get_sequence(project)
        self.project_dna = dict()
        self.project_dna_sax = dict()

    def start_index(self):
        return 0

    def length(self):
        return len(self.sequence)

    def index_from_int_to_date(self, index):
        return self.sequence.index[index]

    def index_date_to_int(self, index_date):
        return pd.to_datetime(self.sequence.index).get_loc(index_date).start

    def _create_sax_signature(self, sax):
        return "%d_%d_%d" % (sax.window, sax.nbins, sax.nlevels)

    def transform_project_dna(self, sax):
        tranformed_seq = []
        signature = self._create_sax_signature(sax)
        self.project_dna_sax[signature] = sax
        window_index = sax.sliding_window_index(len(self.sequence))
        for idx in window_index:
            start_index = idx.start
            end_index = idx.stop - 1
            window = self.sequence[idx]
            word = sax.symbolize_window(window)
            seq_window = ProjectDNASequenceWindow(self, self.sequence,
                                                  word, start_index, end_index, sax)
            tranformed_seq.append(seq_window)

        self.project_dna[signature] = tranformed_seq


class AlignmentMotif(object):

    def __init__(self):
        self.alignments = set()
        self.words = set()

    def word_belongs_to_motif(self, word):
        return word in self.words

    def alignment_belongs_to_motif(self, alignment):
        return alignment in self.alignments

    def add_alignment(self, alignment):
        if not self.alignment_belongs_to_motif(alignment):
            self.words.add(alignment.window_a.word)
            self.words.add(alignment.window_b.word)
            self.alignments.add(alignment)

    def __repr__(self):
        return "|".join(sorted(list(self.words)))

    def __str__(self):
        return self.__repr__()

    def __unicode__(self):
        return self.__repr__()


class Alignment(object):
    def __init__(self, window_a, window_b, alignment, similarity_perc_score):
        self.window_a = window_a
        self.window_b = window_b
        self.alignment = alignment
        self.alignment_score = alignment[2]
        self.similarity_perc_score = similarity_perc_score
        self.alignment_start_index = alignment[3]
        self.alignment_end_index = alignment[4]

    def projects(self):
        return (self.window_a.project, self.window_b.project)

    def projects_seq(self):
        return (self.window_a.project_sequence, self.window_b.project_sequence)

    def __repr__(self):
        return str((self.alignment, self.similarity_perc_score))

    def __eq__(self, other):
        return (isinstance(other, type(self)) and
                (self.window_a, self.window_b) ==
                (other.window_a, other.window_b))

    def __hash__(self):
        return hash((self.window_a, self.window_b))


#############




class ProjectDNASequenceWindow(object):
    def __init__(self, project_sequence, sequence, word, start_index, end_index, sax):
        self.project = project_sequence.project
        self.project_sequence = project_sequence
        self.sequence = sequence
        self.start_index = start_index
        self.end_index = end_index
        self.word = word
        self.window = sax.window
        self.paa = sax.nbins
        self.alphabet = sax.alphabet
        self.alphabet_size = sax.nlevels
        self.sax = sax
        #self.score_func, self.gap_func = create_score_functions(self.sax)


    def __repr__(self):
        return str((self.word, self.start_index, self.end_index))

    def __eq__(self, other):
        return (isinstance(other, type(self)) and
                (self.word, self.window, self.paa, self.alphabet, self.start_index, self.end_index) ==
                (other.word, other.window, other.paa, other.alphabet, self.start_index, self.end_index))

    def __hash__(self):
        return hash((self.word, self.window, self.paa, "".join(self.alphabet), self.start_index, self.end_index))


class ProjectSequence(object):
    def __init__(self, project, get_sequence):
        self.project = project
        self.sequence = get_sequence(project)
        self.project_dna = dict()
        self.project_dna_sax = dict()

    def start_index(self):
        return 0

    def length(self):
        return len(self.sequence)

    def index_from_int_to_date(self, index):
        return self.sequence.index[index]

    def index_date_to_int(self, index_date):
        return pd.to_datetime(self.sequence.index).get_loc(index_date).start

    def _create_sax_signature(self, sax):
        return "%d_%d_%d" % (sax.window, sax.nbins, sax.nlevels)

    def transform_project_dna(self, sax):
        tranformed_seq = []
        signature = self._create_sax_signature(sax)
        self.project_dna_sax[signature] = sax
        # windows = sax.symbolize_signal(self.sequence)
        window_index = sax.sliding_window_index(len(self.sequence))
        for idx in window_index:
            start_index = idx.start
            end_index = idx.stop - 1
            window = self.sequence[idx]
            word = sax.symbolize_window(window)
            seq_window = ProjectDNASequenceWindow(self, self.sequence,
                                                  word, start_index, end_index, sax)
            tranformed_seq.append(seq_window)

        self.project_dna[signature] = tranformed_seq


class ProjecSubsequence(ProjectSequence):

    def __init__(self, project_sequence, start_index, end_index):
        self.project = project_sequence.project
        self.sequence = project_sequence.sequence[start_index:end_index]
        self.project_dna = project_sequence.project_dna
        self.project_dna_sax = project_sequence.project_dna_sax
        self.in_parent_start = start_index
        self.in_parent_end = end_index
        self.parent = project_sequence

    def start_index(self):
        return self.in_parent_start

