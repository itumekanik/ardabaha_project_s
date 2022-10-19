#coding:utf8
import time

_import_time_ = time.time()
_AFT_ = 'TOTAL'

class TIMER_META(type):
    def __getitem__(cls, item):
        if len(cls.clocks)==0: cls.register(_AFT_)
        clock = cls.clocks.get(item, None)
        if clock is None: return cls.register(item)
        if clock.running: clock.stop
        else: clock.start

class clock:
    def __init__(self, name):
        self.name = name
        self.data = [time.time(), 0.0, True] #[start, dT, on/off]

    @property
    def stop(self):
        self.data[2] = False
        self.data[1] += time.time() - self.data[0];

    @property
    def start(self):self.data = [time.time(), self.data[1], True]

    @property
    def checkout(self):
        if self.running:
            self.data[1] += time.time() - self.data[0];
            self.data[0] = time.time()

    @property
    def running(self): return self.data[2]

    @property
    def elapsed(self): return self.data[1]


class timeit(object, metaclass=TIMER_META):
    clocks = {}

    @classmethod
    def register(cls, name): cls.clocks[name] = clock(name)   #[start, dT, on/off]

    @classmethod
    def report(cls, base_clock_name=None):
        print("\n")
        print("-"*10+" Timing Report "+"-"*10)
        for c in cls.clocks.values(): c.checkout
        clocks = sorted(cls.clocks.values(),key=lambda x: x.elapsed, reverse=True)
        if base_clock_name is None: total = cls.clocks[_AFT_].elapsed
        else: total = cls.clocks[base_clock_name].elapsed

        if total==0 and not base_clock_name is None:
            print("Warning: Clock '", base_clock_name, "' elapsed almost zero time!!!")
            print("Base clock is set to TOTAL.")
            total = cls.clocks[_AFT_].elapsed

        if total==0:
            print("Unmeasurable time ratios:")
            print("Setting based elapsed time to 1 sec")
            total = 1.0

        for t in clocks:
            print(t.name + ":", str(t.elapsed)[0:5] + " sec ", "%" + str(t.elapsed / total * 100)[0:4] + "  ")
        print("-" * 6 + " End of Timing Report " + "-" * 7)

    @classmethod
    def print(cls, name):
        for n, t in cls.clocks.items():
            if n.upper()==name.upper():
                print(n.upper() + ":", str(t.elapsed)[0:5])


"""
for j in range(4):

    timeit["dene1"]

    a=0
    for i in range(int(1e6)):
        a = a+i


    timeit["dene1"]

    timeit["dene2"]

    a = 0
    for i in range(int(2e6)):
        a = a + i

    timeit["dene2"]

    timeit.report()

    timeit["denea"]
    a=0
    for i in range(int(1e2)):
        a = a+i
    timeit["denea"]


timeit.report()

"""


