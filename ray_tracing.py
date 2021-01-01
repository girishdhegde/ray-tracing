import numpy as np
import cv2
import time

class point:

    def __init__(self, x=0, y=0, z=0):
        if isinstance(x, np.ndarray):
            self.p = x
        else:
            self.p = np.array([x, y, z])

    def __add__(self, other):
        assert isinstance(other, point)
        return point(self.p+other.p)

    def __sub__(self, other):
        assert isinstance(other, point)
        return point(self.p-other.p)

    def __mul__(self, other):
        # if vector
        if isinstance(other, point):
            return point(self.p*other.p)
        # if scalar
        else:
            return point(self.p*other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
         # if vector
        if isinstance(other, point):
            return point(self.p/other.p)
        else:
            return point(self.p/other)
    
    def __pow__(self, exp):
        return point(self.p ** exp)
    
    def __neg__(self):
        return point(self.p * -1)
        
    def __len__(self):
        return np.linalg.norm(self.p)
    
    def dot(self, other):
        assert isinstance(other, point)
        return self.p.dot(other.p)
    
    def cross(self, other):
        assert isinstance(other, point)
        return point(*np.cross(self.p, other.p))

    def length(self):
        return np.linalg.norm(self.p)
    
    def unitVector(self):
        return point(self.p/np.linalg.norm(self.p))
    
    def getPos(self):
        return (self.p[0], self.p[1], self.p[2])
    
    def distance(self, other):
        if other is None:
            return None
        assert isinstance(other, point)
        return (self - other).length()

    def __str__(self):
        return f'point/vector at {self.p[0], self.p[1], self.p[2]}'
    
    def __repr__(self):
        return f'point/vector at {self.p[0], self.p[1], self.p[2]}'


class vector(point):
    pass


def insideOutside(p, N, *v):
    for i in range(len(v)-1):
        edge = v[i+1] - v[i]
        c = p - v[i]
        if N.dot(edge.cross(c)) < 0:
            return None
    edge = v[0] - v[-1]
    c = p - v[-1]
    if N.dot(edge.cross(c)) < 0:
        return None
    return p


class ray:

    def __init__(self, x=0, y=0, z=0, direction=None):
        self.source = point(x, y, z)
        # Unit vector representing direction
        self.direction = direction

    def __str__(self):
        return f'Ray with origin {self.source} and direction {self.direction}'
    
    def __repr__(self):
        return f'Ray with origin {self.source} and direction {self.direction}'

    def point_at(self, t):
        return self.source + t * self.direction

    def intersects(self, obj):
        if isinstance(obj, sphere):
            x, n = self._raySphereX(obj)
            return x, self.source.distance(x), n
        elif isinstance(obj, triangle):
            x, n = self._rayTriX(obj)
            return x, self.source.distance(x), n
        elif isinstance(obj, rectangle):
            x, n = self._rayRectX(obj)
            return x, self.source.distance(x), n
        elif isinstance(obj, pyramid):
            x, n = self._rayPyrX(obj)
            return x, self.source.distance(x), n
        elif isinstance(obj, cylinder):
            x, n = self._rayCylinderX(obj)
            return x, self.source.distance(x), n
        # elif isinstance(obj, plane):
        #     return self._rayPlaneX(obj)
            
    def _raySphereX(self, obj):
        # a = dot(ray_direction, ray_direction)
        a   = 1
        s_c = self.source - obj.centre
        # b = 2 . dot(ray_direction, sourceToCentre)
        b   = 2 * self.direction.dot(s_c)
        # c = dot(sourceToCentre, sourceToCentre) - radius ** 2
        c   = s_c.dot(s_c) - (radius ** 2)

        discriminant = b ** 2 - 4 * a * c

        # no intersection
        if discriminant < 0:
            return None, None
        # +ve - 2 intersections(passes), 0 - one intersection(tangent)
        else:
            # checking for intersections behind the ray i.e behind camera
            num = -b - discriminant ** 0.5
            if num < 0:
                num = -b + discriminant ** 0.5
                if num < 0:
                    return None, None
            x = self.point_at(num / (2 * a))
            return  x, obj.normal(x)
        
    def _rayPlaneX(self, obj):
        n = obj.normal()
        denom = n.dot(self.direction)
        # normals are opposite in direction
        if denom < -1e-6:
            t = ((obj.p0 - self.source).dot(n)) / denom
            return self.point_at(t), n
        if denom > 1e-6:
            t = ((obj.p0 - self.source).dot(n)) / denom
            return self.point_at(t), n
        return None, None
    
    def _rayRectX(self, obj):
        p, _ = self._rayPlaneX(obj.plane)
        if p is None:
            return None, None
        # inside - outside
        p = insideOutside(p, obj.norm, obj.v0, obj.v1, obj.v2, obj.v3)

        if p is None:
            return None, None
        
        # check if triangle is behind the ray
        toP = p - self.source
        if toP.dot(self.direction) < 0:
            return None, None

        return p, obj.norm     

    def _rayTriX(self, obj):
        p, _ = self._rayPlaneX(obj.plane)
        if p is None:
            return None, None
        # inside - outside
        p = insideOutside(p, obj.norm, obj.v0, obj.v1, obj.v2)

        if p is None:
            return None, None
        
        # check if triangle is behind the ray
        toP = p - self.source
        if toP.dot(self.direction) < 0:
            return None, None

        return p, obj.norm     

    def _rayPyrX(self, obj):
        for face in obj.faces:
            p, n = self._rayTriX(face)
            if p is not None:
                # obj.norm = face.norm
                return p, n
        p, n = self._rayRectX(obj.base)
        # obj.norm = obj.base.norm
        return p, n
    
    def _rayCylinderX(self, obj):
        va = obj.va
        pa = obj.pa
        p1 = pa
        p2 = obj.p2
        r  = obj.r

        # intersection with lateral surface
        delP  = self.source - pa
        temp1 = self.direction - (self.direction.dot(va) * va)
        temp2 = delP - delP.dot(va) * va
        
        a = temp1.dot(temp1)
        b = 2 * temp1.dot(temp2)
        c = temp2.dot(temp2) - (r ** 2)
        
        discr = b ** 2 - 4 * a * c
        if discr < 0:
            t1, t2 =  None, None
        else:
            t1 = (-b + discr ** 0.5) / (2 * a)
            t2 = (-b - discr ** 0.5) / (2 * a)

            if t1 >= 0:
                q1 = self.point_at(t1)
                if va.dot(q1 - p1) > 0 and va.dot(q1 - p2) < 0:
                    pass
                else: 
                    t1 = None
            else:
                t1 = None

            if t2 >= 0:
                q2 = self.point_at(t2)
                if va.dot(q2 - p1) > 0 and va.dot(q2 - p2) < 0:
                    pass
                else: 
                    t2 = None
            else:
                t2 = None
        t1 = None
        # intersection with bottom disc
        t3 = None
        # x, n = self._rayPlaneX(obj.bottomPlane)
        # if x is not None:
        #     xp = x - p1
        #     if xp.dot(xp) < (r ** 2):
        #         t3 = (x - self.source) / self.direction
        #         n3 = n 

        # intersection with top disc
        t4 = None
        # x, n = self._rayPlaneX(obj.topPlane)
        # if x is not None:
        #     xp = x - p2
        #     if xp.dot(xp) < (r ** 2):
        #         t4 = (x - self.source) / self.direction
        #         n4 = n

        t = min([t1 if t1 is not None else float('inf'), t2 if t2 is not None else float('inf'), 
                 t3 if t3 is not None else float('inf'), t4 if t4 is not None else float('inf')])
        x = None
        n = None
        if t < float('inf'):
            x = self.point_at(t) 
            n = obj.normal(x)
        
        return x, n


class light:

    def __init__(self, x=0, y=0, z=0, clr=(255, 255, 255)):
        self.position = point(x, y, z)
        self.clr = np.array(clr)


class material:

    def __init__(self, clr=(255, 255, 255), ambient=0.1, diffuse=0.5, specular=0.5, 
                 reflection=0.3, k=32, ri=0.0, transparency=None):
        self.clr = np.array(clr)
        self.ambient  = ambient
        self.diffuse  = diffuse
        self.specular = specular
        self.reflection = reflection
        self.transparency = transparency
        self.ri = ri
        self.k  = k
    
    def getColor(self, position):
        return self.clr


class floorMaterial:

    def __init__(self, clr1=(255, 255, 255), clr2=(1, 1, 1), tileSize=100, ambient=0.1,
                 diffuse=0.6, specular=1.0, reflection=0.9, k=32, ri=0, transparency=None):
        self.clr1     = np.array(clr1)
        self.clr2     = np.array(clr2)
        self.ambient  = ambient
        self.diffuse  = diffuse
        self.specular = specular
        self.tileSize = tileSize
        self.reflection = reflection
        self.transparency = transparency
        self.ri = ri
        self.k  = k

    def getColor(self, position):
        x, _, z = position.p
        # if odd, odd or even, even
        if int((x // self.tileSize) % 2) ^ int((z // self.tileSize) % 2):
            return self.clr1
        return self.clr2


class sphere:
    def __init__(self, centre, radius, material=None):
        self.centre = point(*centre)
        self.radius = radius
        self.material = material

    def normal(self, position):
        return (position - self.centre).unitVector()


class plane:
    def __init__(self, pt0, pt1=None, pt2=None, norm=None):
        self.p0  = point(*pt0)
        if norm is not None:
            if isinstance(norm, point):
                self.norm = norm
            else:
                self.norm = point(*norm)
        else:
            v0v1 = vector(*pt1) - self.p0
            v0v2 = vector(*pt2) - self.p0
            self.norm = (v0v1.cross(v0v2)).unitVector()

    def normal(self, *args):
        return self.norm


class triangle:
    def __init__(self, v0, v1, v2, material=None):
        '''
        v0, v1, v2 are vertices of triangle taken in CCW order
        ''' 
        self.v0 = point(*v0)
        self.v1 = point(*v1)
        self.v2 = point(*v2)
        self.material = material
        self.plane = plane(v0, v1, v2)
        self.norm  = self.plane.norm

    def normal(self, *args):
        return self.norm


class rectangle:
    def __init__(self, v0, v1, v2, v3, material=None):
        '''
        v0, v1, v2, v3 are vertices of rectangle 
        '''   
        self.v0 = point(*v0)
        self.v1 = point(*v1)
        self.v2 = point(*v2)
        self.v3 = point(*v3)
        self.material = material
        self.plane = plane(v0, v1, v2)
        self.norm  = self.plane.norm

    def normal(self, *args):
        return self.norm


class pyramid:
    def __init__(self, v0, v1, v2, material=None):
        '''
        faceTriangle: leftVertex, rightVertex, topVertex
        '''
        faceTriangle = triangle(v0, v1, v2, material)
        # perpendiculars
        norm = faceTriangle.norm
        p0 = point(*v0)
        p1 = point(*v1)
        dist = p0.distance(p1)
        tmp0 = (p0 + dist * norm).p
        tmp1 = (p1 + dist * norm).p
        self.base  = rectangle(  v0,   v1, tmp1, tmp0, material)
        self.faces = [triangle(  v0,   v1, v2, material),
                      triangle(  v1, tmp1, v2, material),
                      triangle(tmp1, tmp0, v2, material),
                      triangle(tmp0,   v0, v2, material)]
        self.material = material
        self.norm = None
    
    def normal(self, pos):
        return self.norm



class cylinder:
    def __init__(self, p1, p2, r, material=None):
        '''
        p1 - position vector describing the first end point of the long axis of the cylinder 
        p2 - position vector describing the second end point of the long axis of the cylinder 
        r  - radius 
        '''
        self.p1 = point(*p1)
        self.p2 = point(*p2)
        self.r  = r
        self.pa = self.p1
        # unit vector along the axis of the cylinder
        self.va = (self.p2 - self.p1).unitVector()
        # bottom disc
        self.bottomPlane = plane(p1, norm=-self.va)
        # top disc
        self.topPlane = plane(p2, norm=self.va)
        self.material = material
    
    def normal(self, pos):
        # if point in bottom disc
        # temp = pos - self.p1
        # if temp.dot(temp) < (self.r ** 2):
        #     return -self.va
        # # if point in top disc
        # temp = pos - self.p2
        # if temp.dot(temp) < (self.r ** 2):
        #     return self.va
        # if point in lateral surface
        return  self.va.dot(pos) * self.va
        

class box:
    pass

class cone:
    pass


class scene:
    def __init__(self, camera=(0, 0, -100), objs=None, lights=None):
        self.objs   = objs if objs is not None else []
        self.camera = point(*camera)
        self.lights = lights

    def addObject(self, obj):
        self.objs.append(obj)
    
    def addLight(self, light):
        self.lights.append(light)


class engine:
    '''
    function to get image(open-cv pixel array) co-ordinates from worldco-ordinates
    ''' 
    def _getPos(self, x, y, w, h):
        return (h - (y + h // 2), x + w // 2)

    def render(self, scene, width, height, bg=(0, 0, 0), maxRec=5, show=True, savePath=None, shadow=True):
        start = time.time()
        self.maxRec = maxRec
        self.shadow = shadow
        screen = np.zeros((height, width, 3), dtype=np.uint8)
        print('RENDERING IN PROGREESS')
        for x in range(-width//2+1, width//2):
            for y in range(-height//2+1, height//2):
                pos = self._getPos(x, y, width, height)
                pt  = point(x, y, 0)
                r   = ray(scene.camera.p, direction=(pt - scene.camera).unitVector())
                clr = self.rayTrace(r, scene)
                if clr is not None:
                    screen[pos] = tuple(clr)
            over = int(1+50*(2*x+width)//width)
            dash = '-' * (over - 1)
            star = ' ' * (100 - over)
            print('[', dash, '>', star, ']', over, '%', end='\r')
        print('\nrendered in', time.time()-start, 'sec')
        if show:
            cv2.imshow('rayTracing', screen)
            cv2.waitKey(0)
            cv2.destroyAllWindows()
        if savePath is not None:
            cv2.imwrite(savePath, screen)

    def rayTrace(self, r, scene, rec=1):
        if rec > self.maxRec:
            return None
        camera = scene.camera
        minObj, minDist, hit, normal = self.nearest(scene.objs, r)
        if minObj is not None:
            if normal.dot(r.direction) > 0 and not (isinstance(minObj, sphere)):
                normal = -normal
            cameraVec = (camera - hit).unitVector()
            # direct componet
            dirClr = self.phongShader(scene, minObj.material, normal, scene.lights, hit, cameraVec)
            # if refractive
            ri = minObj.material.ri
            tr = minObj.material.transparency
            if ri > 0:
                refractionColor = 0; 
                # compute fresnel
                kr, kt = self.fresnel(r.direction, normal, ri) 
                #############33
                outside = r.direction.dot(normal) < 0 
                bias    = normal * 1e-4
                # refraction if no  TIR
                if kr < 1:
                    refractionDirection = self.refract(r.direction, normal, ri)
                    if refractionDirection is not None:
                        refrdir = refractionDirection.unitVector() 
                        refrorg = hit - bias if outside else hit + bias 
                        refrclr = self.rayTrace(ray(refrorg.p, direction=refrdir), scene, rec+1)
                        if refrclr is not None:
                            refrClr = refrclr * kt * tr
                            dirClr = np.clip(dirClr + refrClr, 0, 255) 
            
                # Reflection
                refPos = hit + normal * 1e-4
                refDir = (r.direction - 2 * (r.direction.dot(normal)) * normal).unitVector()
                refRay = ray(refPos.p, direction=refDir)
                reflc  = self.rayTrace(refRay, scene, rec+1)
                if reflc is not None:
                    refClr = reflc * kr * tr
                    dirClr = np.clip(dirClr + refClr, 0, 255) 
            elif minObj.material.reflection > 0:
                # Reflection
                refPos = hit + normal * 1e-4
                refDir = (r.direction - 2 * (r.direction.dot(normal)) * normal).unitVector()
                refRay = ray(refPos.p, direction=refDir)
                reflc  = self.rayTrace(refRay, scene, rec+1)
                if reflc is not None:
                    refClr = minObj.material.reflection * reflc
                    dirClr = np.clip(dirClr + refClr, 0, 255)   
            return  dirClr
        return None
    
    def nearest(self, objs, ray):
        minDist = float('inf')
        minObj  = None
        hit = None
        normal = None
        for obj in objs:
            intr, dist, norm = ray.intersects(obj)
            if (intr is not None) and (dist < minDist):
                minDist = dist
                minObj  = obj
                hit     = intr
                normal  = norm
        return minObj, minDist, hit, normal

    def phongShader(self, scene, material, normal, lights, hit, cameraVec):
        mat    = material
        amb    = mat.ambient
        dif    = mat.diffuse
        spc    = mat.specular
        spk    = mat.k
        clr    = mat.getColor(hit)
        ambClr = np.array([0, 0, 0])
        ambClr = ambClr + amb * clr
        difClr = np.array([0, 0, 0])
        spcClr = np.array([0, 0, 0])
        for L in lights:
            lightVec = (L.position - hit).unitVector()  
            # shadow
            if self.shadow:
                shifted  = hit + 1e-4 * normal
                minObj, minDist, _, _ = self.nearest(scene.objs, ray(shifted.p, direction=lightVec))
                if minDist is not None and minDist < hit.distance(L.position):
                    # just a hack to represent light through transparency materials
                    if minObj.material.transparency is None:
                        continue
                    else:
                        continue
                        spc *=  (minObj.material.transparency / 4)
                        dif *=  (minObj.material.transparency / 4)
                
            # ambClr = ambClr + amb * clr
            # Lambert shading - diffuse reflection
            nL = max(normal.dot(lightVec), 0)
            difClr = np.clip(difClr+dif*clr*nL, 0, 255) 
            # Blin-Phong shading - specular reflection
            half = (lightVec + cameraVec).unitVector()
            nH   = (max(normal.dot(half), 0)) ** spk
            spcClr = np.clip(spcClr+spc*L.clr*nH, 0, 255)

        return np.clip(ambClr+difClr+spcClr, 0, 255)

    def fresnel(self, i, n, ri):
        ''' 
        i  - unit incident ray dircetion vector
        n  - unit normal vector to the plane
        ri - refractive index of the medium
        '''       
        cosi = np.clip(i.dot(n), -1, 1)
        etai = 1
        etat = ri 
        if cosi > 0:
            etai, etat = etat, etai 
        # Snell's law
        sint = (etai / etat) * (max(0, 1 - cosi * cosi)) ** 0.5 
        # TIR
        if sint >= 1: 
            kr = 1
        else:
            cost = (max(0, 1 - sint * sint)) ** 0.5 
            cosi = abs(cosi) 
            Rs   = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost)) 
            Rp   = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost)) 
            kr   = (Rs * Rs + Rp * Rp) / 2; 

        return kr, 1 - kr


    def refract(self, i, n, ri): 
        ''' 
        i  - unit incident ray dircetion vector
        n  - unit normal vector to the plane
        ri - refractive index of the medium
        '''       
        cosi = np.clip(i.dot(n), -1, 1)
        etai = 1
        etat = ri 
        norm = n
        if cosi < 0:
            cosi = -cosi
        else: 
            etai, etat = etat, etai
            norm = -n 
        eta = etai / etat 
        k = 1 - eta * eta * (1 - cosi * cosi)
        if k < 0:
            return None 
        return eta * i + (eta * cosi - (k ** 0.5)) * norm 
 
if __name__ == '__main__':
    light1 = light(200, 400, 10, (255, 255, 255))
    # red reflective ball
    ball1  = sphere((-100,    0, 500), 200, material((  0,   0, 255), 0.1, 0.8, 1.0, 0.5, 200))
    # purple ball
    ball2  = sphere((-600, -100, 650), 100, material((255,   0, 150), 0.1, 0.6, 0.4, 0.15, 16))
    # yellow semi transparency ball
    ball3  = sphere(( 250,  -50, 200), 150, material(( 50, 255, 255), 0.3, 0.0, 0.0, 0.0,  0, 1.5, 0.5))
   
    floor = sphere((0, -40000, 0), 40000-200, floorMaterial(reflection=0.3))
    scn   = scene((0, 0, -400), [ball1, ball2, ball3, floor, ], [light1, ])
    eng   = engine()
    eng.render(scn, 600, 400, (255, 255, 255), savePath='./balls.png', shadow=True, maxRec=5)

    # tri0  = triangle((-550, -200, 400), (-400, -200, 600), (-475, 200, 500), material((0, 50, 150), 0.1, 0.8, 0.8, 0.5, 32))

    # mirror
    # rect0 = rectangle((-550, -200, 700), (550, -200, 700), (550, 400, 700), (-550, 400, 700),
    #                   material((192, 192, 192), 0.1, 0.5, 0.8, 0.9, 32))
    
    # glass ball in room

    # glass = material((198, 226, 227), ambient=0.01, diffuse=0.0, specular=1.0, 
    #                   reflection=0.0, k=200, ri=1.5, transparency=.9)
    # rect0  = sphere(( 0,  -100, 200), 100, material=glass)

    # floor = sphere((0, -40000, 0), 40000-600, floorMaterial(reflection=0.3))
    # # back  = sphere((0, 0, 80000), 60000, material((255, 255, 255), 0.1, 0.9, 0.0, 0.0, 0))
    # wall1 = sphere((-40000, 0, 0), 40000-800, material((0, 255, 255), 0.01, 0.9, 0.0, 0.0, 0))
    # wall2 = sphere((40000, 0, 0), 40000-800, material((0, 0, 255), 0.01, 0.9, 0.0, 0.0, 0))
    # ceil  = sphere((0, 40000, 0), 40000-600, material((255, 255, 255), 0.3, 0.9, 0.0, 0.0, 0))
    # scn   = scene((0, 50, -400), [floor, wall1, wall2, ceil, rect0  ], [light1, ])
    
    # To be updated
    # pyr0  = pyramid((-500, -200, 400), (0, -200, 200), (-250, 200, 300), material((0, 50, 150), 0.1, 0.9, 0.1, 0.1, 32))
    # scn   = scene((0, 0, -400), [pyr0, floor, ball], [l1, ])
    # cyl   = cylinder((0, -200, 400), (0, 200, 400), 100, material((0, 0, 255), 0.2, 0.3, 0.0, 0.0, 0))
    # scn   = scene((0, 50, -400), [cyl, ], [light1, ])
    # eng   = engine()
    # eng.render(scn, 600, 400, (255, 255, 255), savePath='./balls.png', shadow=True, maxRec=5)
    # l2    = light(200, 200, 10, (255, 255, 255))

    
    '''
    To Do:
    *. movable screen
    *. soft shadow
    *. glass

    *. other shapes: cone, cylinder, cuboid
    *. Textures
    *. mesh rendering from .obj files
    '''
