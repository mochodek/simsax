import pandas as pd
import os

class Project(object):
    def __init__(self, name, source,  dev_type, backlog):
        self.name = name
        self.source = source
        self.dev_type = dev_type
        self.backlog = backlog

    def from_date(self):
        return self.backlog.from_date.min()

    def to_date(self):
        return self.backlog.from_date.max()

    def id(self):
        return "%s-%s" % (self.source, self.name)

    def __unicode__(self):
        return "<%s, %s, %s>" % (self.name, self.source, self.dev_type)

    def __str__(self):
        return self.__unicode__()

    def __repr__(self):
        return self.__unicode__()


def load_backlog(backlogs_folder, file_name, source, dev_type):
    file_path = backlogs_folder+file_name
    defect_backlog = pd.read_csv(file_path, parse_dates=['from_date', 'to_date'], index_col=0)
    name = file_name.split('-')[0]
    project = Project(name=name, source=source, dev_type=dev_type, backlog=defect_backlog)
    return project

