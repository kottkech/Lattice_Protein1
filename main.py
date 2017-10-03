import time

class Protein:
    def generateRes(self, angles):
        for angle in angles:
            len = self.length
            r1 = self.chain[len-2]
            r2 = self.chain[len-1]
            dx = r2.x - r1.x
            dy = r2.y - r1.y

            nx = 0
            ny = 0
            if angle == 0:
                nx = dx
                ny = dy
            elif angle == -1:
                if dx == 0:
                    ny = 0
                    nx = -1 * dy
                if dy == 0:
                    nx = 0
                    ny = dx
            elif angle == 1:
                if dx == 0:
                    ny = 0
                    nx = dy
                if dy == 0:
                    nx = 0
                    ny = -1 * dx

            nx = nx + r2.x
            ny = ny + r2.y

            if self.check(nx,ny):
                self.chain.append(Res(nx,ny))
                self.length += 1
            else:
                raise ValueError("The given protein self-intersects")

    def check(self, x, y):

        for res in self.chain:
            if res.x == x and res.y == y:
                return False

        return True

    def __init__(self, angles):
        self.chain = [Res(0,0), Res(0,1)]
        self.length = 2
        self.generateRes(angles)

class Proteins:
    def gyration(self):
        aveGy = 0
        minGy = -1
       
        for chain in self.chains:
            straight = False
            mx = 0
            my = 0
            sum = 0
            for res in chain:
                mx += res.x
                my += res.y
                if(res.y == len(chain)-1):
                    straight = True
                    
            mx /= float(len(chain))
            my /= float(len(chain))

            for res in chain:
                sum += ((res.x - mx)**2 + (res.y - my)**2)

            gy = sum / float(len(chain))
            gy = gy**(1/2)

            if straight:
                aveGy += (4*gy) #to account for the 4-fold rotational symmetry
            else:
                aveGy += (8*gy) #to account for the 4-fold rotational symmetry and reflection
                
            if minGy == -1 or gy < minGy:
                minGy = gy

        aveGy /= float(8*len(self.chains) - 4) #accounting for 4 rotations, reflection minus 4 for the non-reflection of the straight chains

        print("Min Gyration: " + str(minGy))
        print("Ave Gyration: " + str(aveGy))


    def add(self, dir, chain):
        size = len(chain)
        r1 = chain[size - 2]
        r2 = chain[size - 1]
        dx = r2.x - r1.x
        dy = r2.y - r1.y

        nx = 0
        ny = 0
        if dir == 0:
            nx = dx
            ny = dy
        elif dir == -1:
            if dx == 0:
                ny = 0
                nx = -1 * dy
            if dy == 0:
                nx = 0
                ny = dx
        elif dir == 1:
            if dx == 0:
                ny = 0
                nx = dy
            if dy == 0:
                nx = 0
                ny = -1 * dx

        nx = nx + r2.x
        ny = ny + r2.y

        return chain + [Res(nx, ny)]

    def intersects(self, chain):
        last = chain[len(chain)-1]
        for i in range(0,len(chain)-1):
            if last.x == chain[i].x and last.y == chain[i].y:
                return True
        return False

    def recGen(self, N, chain, right):
        if N == 0:
            return [chain]
        else:
            out = []

            if right:
                b1 = self.add(-1,chain)
                if not self.intersects(b1):
                    out += self.recGen(N-1,b1, right)
            b2 = self.add(0,chain)
            if not self.intersects(b2):
                out += self.recGen(N-1,b2, right)
            b3 = self.add(1,chain)
            if not self.intersects(b3):
                out += self.recGen(N-1,b3, True)
                

            return out

    def __init__(self, size):
        chain = [Res(0,0), Res(0,1)]
        self.chains = self.recGen(size-2, chain, False)


class Res:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def toString(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"


#p = Protein([1,1,0,0,1,-1])

for i in range(4,17):
    print("i=" + str(i))
    start = time.clock()*1000000
    p = Proteins(i)
    mid = time.clock()*1000000
    p.gyration()
    end = time.clock()*1000000

    print("Total non-intersecting chains: " + str(8*len(p.chains)-4))
    print("Generation: " + str(mid-start) + " micros")
    print("Gyration  : " + str(end-mid) + " micros")
    print("Total     : " + str(end-start) + " micros")
    print("-----")

#3.
# this code has the complexity t=0.6855e^(0.9806N) in microseconds (R^2 = 0.998)
# this means that it would take 2,354.9 years to compute N = 40
# If N = 100, it would take 8.397 * 10^28 years to compute
# The significance of this is that this extremely simple model cannot handle relatively small values of N.
# The one optimization that could be made is if I used a lookup table grid instead of a chain of polymers as this would make looking up individual units O(1) instead of O(n)
# One could speed up the search by pruning chains with three consecutive same-handed turns. This would speed up the process by 22.2% (2/9) as there are 2 such handed turns every three steps.

#4.
#Plotting <R_g> against N, I get the fit <R_g> = 0.3119x^0.7614 . This exponent value is very close to the expected value of 0.75. The difference between the two values could be due to the relatively small sample size.
#Given this value, the lattice model seems to be useful for investigating <R_g>

#Plotting R_min against N, I get the fit R_min = 0.3563x^0.5412. This is near what I would expect as the minimal R value occurs when the chain coils up on itself thus the R_g would be constant. Since side length is N^(1/2) This would mean that the R value would be around N^(1/2)
#Estimating, R_min for N=25 = 1.78 and R_min for N=100 = 3.56.
#After calculating the integral N=25 should = 2.04 and N=100 should = 4.08
