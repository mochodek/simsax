class AlignmentsCoverageCalculator(object):
    def __init__(self, between_alignments, within_alignments=None):
        self.between_alignments = between_alignments
        self.within_alignments = within_alignments
        if len(self.between_alignments) > 0:
            self.projects = self.between_alignments[0].projects()
            self.projects_seq = self.between_alignments[0].projects_seq()

            self.between_covered = {self.projects[0]: list(), self.projects[1]: list()}
            self.between_covered_set = {self.projects[0]: set(), self.projects[1]: set()}

            self.within_covered = {self.projects[0]: list(), self.projects[1]: list()}
            self.within_covered_set = {self.projects[0]: set(), self.projects[1]: set()}

            self.metrics = {self.projects[0]: dict(), self.projects[1]: dict()}
        else:
            self.projects = []
            self.projects_seq = []
            self.metrics = {}
            self.between_covered = {}
            self.between_covered_set = {}
            self.within_covered = {}
            self.within_covered_set = {}


    def calculate_weeks(self):
        self.calculate_between_weeks()
        self.calculate_within_weeks()

    def calculate_between_weeks(self):

        for alignment in self.between_alignments:
            # calculate indexes covered by alignments
            cov_a = [x for x in range(alignment.window_a.start_index, alignment.window_a.end_index + 1)]
            cov_b = [x for x in range(alignment.window_b.start_index, alignment.window_b.end_index + 1)]
            self.between_covered[alignment.window_a.project].extend(cov_a)
            self.between_covered[alignment.window_b.project].extend(cov_b)

        # calculate unique weeks covered
        for project in self.projects:
            self.between_covered_set[project] = set(self.between_covered[project])

        self._calculate_metrics()

    def calculate_within_weeks(self):
        for project in self.projects:
            for alignment in self.within_alignments[project]:
                # calculate indexes covered by alignments
                cov_a = [x for x in range(alignment.window_a.start_index, alignment.window_a.end_index + 1)]
                cov_b = [x for x in range(alignment.window_b.start_index, alignment.window_b.end_index + 1)]
                self.within_covered[project].extend(cov_a)
                self.within_covered[project].extend(cov_b)

            # calculate unique weeks covered
            self.within_covered_set[project] = set(self.within_covered[project])

        self._calculate_metrics()


    def get_novelty_within_weeks(self, project):
        return self.within_covered_set[project] - self.between_covered_set[project]

    def get_novelty_between_weeks(self, project):
        return self.between_covered_set[project] - self.within_covered_set[project]


    def _calculate_metrics(self):
        # calculate common coverage
        for project_seq in self.projects_seq:
            project = project_seq.project
            self.metrics[project]['between_coverage'] = len(self.between_covered_set[project]) / len(
                project_seq.sequence)
            self.metrics[project]['within_coverage'] = len(self.within_covered_set[project]) / len(project_seq.sequence)
            self.metrics[project]['unique_within_coverage'] = len(
                self.get_novelty_within_weeks(project)) / len(project_seq.sequence)
            self.metrics[project]['unique_between_coverage'] = len(
                self.get_novelty_between_weeks(project)) / len(
                project_seq.sequence)



