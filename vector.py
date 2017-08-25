import math

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def plus(self, v):
        new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def minus(self, v):
        new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def times_scalar(self, c):
        new_coordinates = [x*c for x in self.coordinates]
        return Vector(new_coordinates)


    def magnitude(self):
        sumOfSqrs = [x*x for x in self.coordinates]
        return math.sqrt(sum(sumOfSqrs))


    def direction(self):
        mag = self.magnitude()
        if mag == 0:
            print("No direction")
            return None
        new_coordinates = [x/mag for x in self.coordinates]
        return Vector(new_coordinates)


    def dot_product(self, v):
        products = [x*y for x, y in zip(self.coordinates, v.coordinates)]
        return sum(products)


    def angle(self, v):
        mag1 = self.magnitude()
        mag2 = v.magnitude()
        if (mag1==0) or (mag2==0):
            print("No angle")
            return None
        val = self.dot_product(v)/(mag1*mag2)
        if val>1:
            val = 1
        elif val<-1:
            val = -1
        theta = math.acos(val)
        return  theta


    def parallel(self, v):
        mag1 = self.magnitude()
        mag2 = v.magnitude()
        if (mag1==0) or (mag2==0):
            return True
        theta = self.angle(v)
        theta1 = abs(theta)
        theta2 = abs(theta-math.pi)
        if (theta1<0.0001) or (theta2<0.0001):
            return True
        return False


    def orthognal(self, v):
        mag1 = self.magnitude()
        mag2 = v.magnitude()
        if (mag1==0) or (mag2==0):
            return True
        theta = self.angle(v)
        theta1 = abs(theta-(math.pi/2))
        theta2 = abs(theta-3*(math.pi/2))
        if (theta1<0.0001) or (theta2<0.0001):
            return True
        return False


    def parallel_projection_on(self, v):
        dir = v.direction()
        mag = self.dot_product(dir)
        return dir.times_scalar(mag)


    def orthognal_projection_on(self, v):
        parallel = self.parallel_projection_on(v)
        return self.minus(parallel)


    def cross_product(self, v):
        if (self.dimension == 2):
            self.coordinates = self.coordinates +(0,)
            v.coordinates = v.coordinates +(0,)
            self.dimension += 1
            v.dimension += 1
        x = self.coordinates[1]*v.coordinates[2] - v.coordinates[1]*self.coordinates[2]
        y = -(self.coordinates[0]*v.coordinates[2] - v.coordinates[0]*self.coordinates[2])
        z = self.coordinates[0]*v.coordinates[1] - v.coordinates[0]*self.coordinates[1]
        return Vector([x, y, z])


    def area_of_parallelogram_spanned_with(self, v):
        cross = self.cross_product(v)
        return cross.magnitude()


    def area_of_triangle_spanned_with(self, v):
        cross = self.cross_product(v)
        return (cross.magnitude()/2)


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates

