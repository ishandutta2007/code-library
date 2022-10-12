#nlogn
class SegmentTree:
    def __init__(self, xs):
        self.cnts = defaultdict(int)
        self.total = defaultdict(int)
        self.xs = xs

    def update(self, v, tl, tr, l, r, h):
        if l > r: return
        if l == tl and r == tr:
            self.cnts[v] += h
        else:
            tm = (tl + tr)//2
            self.update(v*2, tl, tm, l, min(r, tm), h)
            self.update(v*2+1, tm+1, tr, max(l, tm+1), r, h)
          
        if self.cnts[v] > 0:
            self.total[v] = self.xs[tr + 1] - self.xs[tl]
        else:
            self.total[v] = self.total[v*2] + self.total[v*2+1]
        return self.total[v]
    
class Solution:
    def rectangleArea(self, rectangles):
        xs = sorted(set([x for x1, y1, x2, y2 in rectangles for x in [x1, x2]]))
        xs_i = {x:i for i, x in enumerate(xs)}

        STree = SegmentTree(xs)
        L = []
        for x1, y1, x2, y2 in rectangles:
            L.append([y1, 1, x1, x2])
            L.append([y2, -1, x1, x2])
        L.sort()

        cur_y = cur_x_sum = area = 0
        
        for y, op_cl, x1, x2 in L:
            area += (y - cur_y) * cur_x_sum
            cur_y = y
            STree.update(1, 0,  len(xs) - 1, xs_i[x1], xs_i[x2]-1, op_cl)
            cur_x_sum = STree.total[1]
            
        return area % (10 ** 9 + 7)