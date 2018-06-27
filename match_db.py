from statistics import mean

class Match_db:
    
    def __init__(self, sn, name, pid, alg_len, e_val, mm, gaps, gaps_o, read):
        self.short_name = sn
        self.name = name
        self.count = 1
        self.pid = [pid]
        self.alg_len = [alg_len]
        self.e_val= [e_val]
        self.missmatch = [mm]
        self.gaps = [gaps]
        self.gap_openings = [gaps_o]
        self.read = [read]
    
    def add_read(self, read_name, pid, alg_len, e_val, mm, gaps, gaps_o):
        self.count += 1
        self.read.append(read_name)
        self.pid.append(pid)
        self.alg_len.append(alg_len)
        self.e_val.append(e_val)
        self.missmatch.append(mm)
        self.gaps.append(gaps)
        self.gap_openings.append(gaps_o)

    def __repr__(self):
        return '''Match_db([{sn},{n}, {c}, {pid}, {al}, {e}, {m}, {g}, {go}, {r}])''' .format(sn=self.short_name, n=self.name, \
                            c=self.count, pid=self.pid, \
                            al=self.alg_len, e=self.e_val, \
                            m=self.missmatch, g=self.gaps, \
                            go=self.gap_openings, r=self.read)

    def __str__(self):
        def preview(lst):
            if len(lst) == 0:
                return ""
            elif len(lst) <= 3:
                return ", ".join(lst)
            else:
                return ", ".join(lst[:3]) + ",..."

        return '''Name: {} \nCount: {} \n%ID-avg: {} \nAlign Len-avg: {} \nE-val-avg: {} \nMissmatch-avg: {} \nGaps-avg: {} \nGap openings-avg: {} \nReads: {}'''.format(self.short_name, self.count, mean(self.pid), mean(self.alg_len), mean(self.e_val), mean(self.missmatch), mean(self.gaps), mean(self.gap_openings), preview(self.read) )
