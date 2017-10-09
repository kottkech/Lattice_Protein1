import time
import random
from scipy.constants import codata
import math

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

            if self.check(nx,ny,self.chain):
                self.chain.append(Res(nx,ny))
                self.length += 1
            else:
                raise ValueError("The given protein self-intersects")

    def check(self, x, y, chain):

        for res in chain:
            if res.x == x and res.y == y:
                return False

        return True

    def __init__(self, angles):
        self.chain = [Res(0,0), Res(0,1)]
        self.length = 2
        self.generateRes(angles)

class RandProtein():
    def getNext(self, dir, chain):
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

        return Res(nx, ny)

    def intersects(self, chain):
        last = chain[len(chain)-1]
        for i in range(0,len(chain)-1):
            if last.x == chain[i].x and last.y == chain[i].y:
                return True
        return False

    def genRand(self, len, chain):
        if len == 0:
            return chain
        else:
            dirs = [-1,0,1]
            random.shuffle(dirs)
            for dir in dirs:
                next = self.getNext(dir, chain)
                if not self.intersects(chain + [next]):
                    res = self.genRand(len-1, chain + [next])
                    if res is not None:
                        return res

            return None

    def toDirs(self):
        res = []
        for i in range(2,len(self.chain)):
            r1 = self.chain[i-2]
            r2 = self.chain[i-1]
            r3 = self.chain[i]
            dx1 = r2.x - r1.x
            dy1 = r2.y - r1.y
            dx2 = r3.x - r2.x
            dy2 = r3.y - r2.y

            if dx1 == dx2 and dy1 == dy2:
                res += [0]
            elif dx1 == 0:
                if dy1 == dx2:
                    res += [1]
                else:
                    res += [-1]
            elif dy1 == 0:
                if dx1 == dy2:
                    res += [-1]
                else:
                    res += [1]

        return res

    def pivot(self):

        dirs = self.toDirs()
        arr = [-1,0,1]
        rand = random.randint(0,len(dirs)-1)
        arr.remove(dirs[rand])
        pick = random.randint(0,1)
        dirs[rand] = arr[pick]

        try:
            p = Protein(dirs)
            self.chain = p.chain
            return True
        except ValueError:
            return False

    def pivotEnergy(self, T):

        dirs = self.toDirs()
        arr = [-1,0,1]
        rand = random.randint(0,len(dirs)-1)
        arr.remove(dirs[rand])
        pick = random.randint(0,1)
        dirs[rand] = arr[pick]

        eOld = self.getEnergy(self.chain)

        try:
            p = Protein(dirs)
            eNew = self.getEnergy(p.chain)
            if T *(eOld - eNew) > math.log(random.random()):
                self.chain = p.chain
                return eNew
            return eOld
        except ValueError:
            return eOld

    def getEnergy(self, chain):
        e = 0
        for i in range(0,len(chain)):
            for j in range(i+1,len(chain)):
                if (chain[i].x - chain[j].x)**2 + (chain[i].y - chain[j].y)**2 == 1:
                    e -= 1
        return e

    def __init__(self, len):
        chain = [Res(0, 0), Res(0, 1)]
        self.length = 2
        self.chain = self.genRand(len-2, chain)

class Proteins:
    def gyration(self):
        aveGy = 0
        minGy = -1
        for chain in self.chains:
            mx = 0
            my = 0
            sum = 0
            for res in chain:
                mx += res.x
                my += res.y
            mx /= float(len(chain))
            my /= float(len(chain))

            for res in chain:
                sum += ((res.x - mx)**2 + (res.y - my)**2)

            gy = sum / float(len(chain))
            gy = gy**(1/2)

            aveGy += gy
            if minGy == -1 or gy < minGy:
                minGy = gy

        aveGy /= float(len(self.chains))

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

    def flip(self, chain):
        result = []
        for res in chain:
            result.append(Res(-1*res.x,res.y))
        return result

    def recGen(self, N, chain):
        if N == 0:
            #if not right:
                return [chain]
            #return [chain] + [self.flip(chain)]
        else:
            out = []
            #if right:
            b1 = self.add(-1,chain)
            if not self.intersects(b1):
                out += self.recGen(N - 1, b1)

            b2 = self.add(0,chain)
            if not self.intersects(b2):
                out += self.recGen(N-1,b2)

            b3 = self.add(1, chain)
            if not self.intersects(b3):
                #right = True
                out += self.recGen(N-1,b3)

            return out

    def __init__(self, size):
        chain = [Res(0, 0), Res(0, 1)]
        self.chains = self.recGen(size - 2, chain)


class Res:
    def __init__(self,x,y):
        self.x = x
        self.y = y

    def toString(self):
        return "(" + str(self.x) + ", " + str(self.y) + ")"



#p = Protein([1,1,0,0,1,-1])

#p1 = Protein([1, -1, 0, -1, -1, -1, 0, -1])

'''
f = open("accRatio.csv", "r+")
f.truncate()
f.close()

f = open("accRatio.csv", "a")
f.write("Size,Accepted\n")
f.close()

sweeps=10

for size in range(10,201):
    acceptance=0

    print(size)

    p = RandProtein(size)
    for i in range(0,size*sweeps):
        if p.pivot():
            acceptance += 1

    f = open("accRatio.csv", "a")
    f.write(str(size) + "," + str(acceptance/(size*sweeps)) + '\n')
    f.close()
    
'''

f = open("energyGo2.csv", "r+")
f.truncate()
f.close()

f = open("energyGo2.csv", "a")
f.write("Beta,<E>\n")
f.close()

sweeps=100
d = 0.01
start=1
end=10

#pTemp = Protein([-1,-1,1,1,-1,-1,0,-1,0,0,1,1,0,0])
pTemp = Protein([0,0,-1,-1,0,0,1,1,0,0,-1,-1,0,0])
p = RandProtein(5)
p.chain = pTemp.chain

beta = start

while beta < end:
    average=0

    print(beta)

    for i in range(0,sweeps*16):
        average += p.pivotEnergy(beta)

    f = open("energyGo2.csv", "a")
    f.write(str(beta) + "," + str(average/(sweeps * 16)) + '\n')
    f.close()

    beta *= (1 + d)

'''
for i in range(4,21):
    shortest = -1
    for j in range(0,10):
        p=None
        start = time.clock()*1000000
        p = Proteins(i)
        end = time.clock()*1000000
        t = end - start
        if shortest == -1 or t < shortest:
            shortest = t

    print("i=" + str(i))
    #print(len(p.chains))
    print(str(shortest) + " nanos")
    #p.gyration()
    print("-----")
'''


#---LAB 1---
#3.
# this code has the complexity t=0.4736e^(1.082x) in microseconds (R^2 = 0.998)
# this means that it would take 93944 years to compute N = 40
# If N = 100, it would take 1.47 * 10^33 years to compute
# The significance of this is that this extremely simple model cannot handle relatively small values of N.
# The one optimization that could be made is if I used a lookup table grid instead of a chain of polymers as this would make looking up individual units O(1) instead of O(n)
# One could speed up the search by pruning chains with three consecutive same-handed turns. This would speed up the process by 22.2% (2/9) as there are 2 such handed turns every three steps.

#4.
#Plotting <R_g> against N, I get the fit <R_g> = 0.3119x^0.7614 . This exponent value is very close to the expected value of 0.75. The difference between the two values could be due to the relatively small sample size.
#Given this value, the lattice model seems to be useful for investigating <R_g>

#Plotting R_min against N, I get the fit R_min = 0.3563x^0.5412. This is near what I would expect as the minimal R value occurs when the chain coils up on itself thus the R_g would be constant. Since side length is N^(1/2) This would mean that the R value would be around N^(1/2)
#Estimating, R_min for N=25 = 1.78 and R_min for N=100 = 3.56.
#After calculating the integral N=25 should = 2.04 and N=100 should = 4.08


#---LAB 2---
'''
1. I note that as chain length increases, the acceptance rate decreases in what seems to be a negative log curve. This would make sense as when chain length increases, the chance that any two residues are overlapping increases.

2. Go Latice
In this case, the minimal energy is -24 corresponding to the starting state of the protein given. Initially, the energy of the system averages around -20.5, but as the temperature is lowered, the energy settles to -23 and remains constant. I attribute the significance of the -23 state to be that this is the energy that is balanced by the greatest entropy.

Random Protein
In the case of the random protein, the initial energy begins higher, and it takes more time for the protein to fold into a stable configuration, but when it settles, it does so at an energy of -24.
When the initial temperature is lowered, the energy of the protein almost instantly drops to -23 and stays there
This is because the temperature is no longer high enough to support a exploration of higher energy states.

3.Go Latice 2
Similar to the first go latice, this run had an initial energy that was high but eventually settled to -24 and stayed there.

'''