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

    def __repr__(self):
        return '''Match_db([{sn},{n}, {c}, {pid}, {al}, {e}, {m}, {g}, {go}, {r}])''' .format(sn=self.short_name, n=self.name, \
                            c=self.count, pid=self.pid, \
                            al=self.alg_len, e=self.e_val, \
                            m=self.missmatch, g=self.gaps, \
                            go=self.gap_openings, r=self.read)

    def __str__(self):
        return "Name: {} \nCount: {}".format(self.short_name, self.count)
