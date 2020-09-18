from abc import abstractmethod

from .parallel_analysis_interface import ProcessorPool


class WorkSplitter:
    def __init__(self, controller, group1, group2):
        self.group1 = group1
        self.group2 = group2
        self.controller = controller

    @classmethod
    def setup(cls, ng1, ng2):
        pp, wg = ProcessorPool.from_sizes(
            [(1, "controller"), (ng1, "group1"), (ng2, "group2")]
        )
        groupc = pp["controller"]
        group1 = pp["group1"]
        group2 = pp["group2"]
        obj = cls(groupc, group1, group2)
        obj.run(wg.name)

    def run(self, name):
        if name == "controller":
            self.run_controller()
        elif name == "group1":
            self.run_group1()
        elif name == "group2":
            self.run_group2()
        else:
            raise NotImplementedError

    @abstractmethod
    def run_controller(self):
        pass

    @abstractmethod
    def run_group1(self):
        pass

    @abstractmethod
    def run_group2(self):
        pass
