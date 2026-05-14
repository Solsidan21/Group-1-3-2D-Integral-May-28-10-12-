from  matplotlib.pyplot import *
import matplotlib.pyplot as plt
import numpy as np

class Mesh:
    def __init__(self,coordinate,element,numbering= None):
        if numbering is None:
            self.coordinate = coordinate 
            self.element = element
        else:
            self.coordinate = coordinate 
            self.element = element
            self.numbering = numbering-1
    def x1(self, i):
        return self.coordinate[self.element[i][0]][0]
    def x2(self, i):
        return self.coordinate[self.element[i][1]][0]
    def x3(self, i):
        return self.coordinate[self.element[i][2]][0]
    def y1(self, i):
        return self.coordinate[self.element[i][0]][1]
    def y2(self, i):
        return self.coordinate[self.element[i][1]][1]
    def y3(self, i):
        return self.coordinate[self.element[i][2]][1]
    def det_J(self):
        return (self.x2(self.numbering)-self.x1(self.numbering))*(self.y3(self.numbering)-self.y1(self.numbering))-(self.y2(self.numbering)-self.y1(self.numbering))*(self.x3(self.numbering)-self.x1(self.numbering))
    def min_angle(self):
        a = [self.x2(self.numbering)-self.x1(self.numbering),self.y2(self.numbering)-self.y1(self.numbering)]
        b = [self.x3(self.numbering)-self.x1(self.numbering),self.y3(self.numbering)-self.y1(self.numbering)]
        c = [self.x3(self.numbering)-self.x2(self.numbering),self.y3(self.numbering)-self.y2(self.numbering)]
        return min(np.arccos((a[0]*b[0]+a[1]*b[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(b[0]**2+b[1]**2))), np.arccos((a[0]*c[0]+a[1]*c[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(c[0]**2+c[1]**2))), np.arccos((c[0]*b[0]+c[1]*b[1])/(np.sqrt(c[0]**2+c[1]**2)*np.sqrt(b[0]**2+b[1]**2))))
    def anropning(self):
        p = []
        for i in range(len(self.element)):
            a = [self.x2(i)-self.x1(i),self.y2(i)-self.y1(i)]
            b = [self.x3(i)-self.x1(i),self.y3(i)-self.y1(i)]
            c = [self.x3(i)-self.x2(i),self.y3(i)-self.y2(i)]
            if min(np.arccos((a[0]*b[0]+a[1]*b[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(b[0]**2+b[1]**2))), np.arccos((a[0]*c[0]+a[1]*c[1])/(np.sqrt(a[0]**2+a[1]**2)*np.sqrt(c[0]**2+c[1]**2))), np.arccos((c[0]*b[0]+c[1]*b[1])/(np.sqrt(c[0]**2+c[1]**2)*np.sqrt(b[0]**2+b[1]**2))))<0.1:
                raise TypeError("Bad Angle")
            else:
                p.append((self.x2(i)-self.x1(i))*(self.y3(i)-self.y1(i))-(self.y2(i)-self.y1(i))*(self.x3(i)-self.x1(i)))
        return p
    def f_dash(self):
        anrop = self.anropning()
        q = []
        def f(x,y):
            return 40*x +20*y
        for i in range(len(anrop)):
            q.append(abs(anrop[i])/6*(f(self.x1(i),self.y1(i))+f(self.x2(i),self.y2(i))+f(self.x3(i),self.y3(i))))
        return sum(q)
    def area(self):
        q = []
        for i in range(len(self.element)):
            q.append(1/2*abs(self.x1(i)*(self.y2(i)-self.y3(i))+self.x2(i)*(self.y3(i)-self.y1(i))+self.x3(i)*(self.y1(i)-self.y2(i))))
        return sum(q)
        

def plot_triangles(folder_path, coordinate_file, element_file):    
    with open(folder_path+"\\"+coordinate_file+".txt", "r") as file:
        lines_list = file.readlines()
        Coord1_x = lines_list[0].split()
        Coord1_y = lines_list[1].split()
        file.close()
    with open(folder_path+"\\"+element_file+".txt", "r") as file:
        lines_list = file.readlines()
        Elementnode1_1 = lines_list[0].split()
        Elementnode1_2 = lines_list[1].split()
        Elementnode1_3 = lines_list[2].split()
        file.close()
    
    
    triangle = [[int(float(a)-1), int(float(b)-1), int(float(c)-1)] for a, b, c in zip(Elementnode1_1, Elementnode1_2, Elementnode1_3)]
    coordinates = [[float(a),float(b)] for a, b in zip(Coord1_x, Coord1_y)]

    
    def x1(i):
        return coordinates[triangle[i][0]][0]
    def x2(i):
        return coordinates[triangle[i][1]][0]
    def x3(i):
        return coordinates[triangle[i][2]][0]
    def y1(i):
        return coordinates[triangle[i][0]][1]
    def y2(i):
        return coordinates[triangle[i][1]][1]
    def y3(i):
        return coordinates[triangle[i][2]][1]
    
    boom = []
    for i in range(len(triangle)): #P′=G+k(P−G) där G = tyngdpunkt och k = scale
        k = 0.9
        Gx = (x1(i)+x2(i)+x3(i))/3
        Gy = (y1(i)+y2(i)+y3(i))/3
        X1 = Gx+k*(x1(i)-Gx)
        X2 =  Gx+k*(x2(i)-Gx)
        X3 =  Gx+k*(x3(i)-Gx)
        Y1 =  Gy+k*(y1(i)-Gy)
        Y2 =  Gy+k*(y2(i)-Gy)
        Y3 =  Gy+k*(y3(i)-Gy)
        Z = [[X1,Y1],[X2,Y2],[X3,Y3]]
        boom.append(Z)
        
    for i in range(len(boom)):
        tri1_x = [boom[i][0][0], boom[i][1][0], boom[i][2][0], boom[i][0][0]]
        tri1_y = [boom[i][0][1], boom[i][1][1], boom[i][2][1], boom[i][0][1]]
        plt.plot(tri1_x, tri1_y)
    
    plt.axis('equal')
    plt.show()
    return coordinates, triangle



folder_path = r"C:\Users\Bruno\Desktop\meshes" # döp om ifall du kör på annan dator

coords = ["coord1", "coord2", "coordinates_dolfin_coarse", "coordinates_unitcircle_400","coordinates_unitcircle_1024", "coordinates_unitcircle_2500","coordinates_unitcircle_10000"]
elements = ["elementnode1", "elementnode2","nodes_dolfin_coarse","nodes_unitcircle_400","nodes_unitcircle_1024", "nodes_unitcircle_2500","nodes_unitcircle_10000"]

index = 0
for lists in zip(coords, elements): #print pairs
    print(index, lists[0], lists[1])
    index += 1

pair = input("What list pair do you want?: ")

coord_list = coords[int(pair)]
element_list = elements[int(pair)]

coordinates, triangles = plot_triangles(folder_path, coord_list, element_list)

v = Mesh(coordinates,triangles)
print(v.f_dash())
print(v.area())