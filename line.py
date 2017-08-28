from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = [0]*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal(0)
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = [0]*self.dimension

            initial_index = Line.first_nonzero_index(n.coordinates)
            initial_coefficient = n.coordinates[initial_index]

            basepoint_coords[initial_index] = float(c)/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def parallel(self, l):
        return self.normal_vector.parallel(l.normal_vector)


    def __eq__(self, l):
        if not(self.parallel(l)):
            return False
        if (self.constant_term==0) and (l.constant_term==0):
            return True
        if (self.constant_term==0) ^ (l.constant_term==0):
            return False
        n1 = self.normal_vector.times_scalar(1/float(self.constant_term))
        n2 = l.normal_vector.times_scalar(1 / float(l.constant_term))
        return (n1 == n2)


    def intersection(self, l):
        if self==l:
            print("Infinite intersections")
            return
        if self.parallel(l):
            print("No intersection")
            return
        A = self.normal_vector.coordinates[0]
        B = self.normal_vector.coordinates[1]
        C = l.normal_vector.coordinates[0]
        D = l.normal_vector.coordinates[1]
        k1 = float(self.constant_term)
        k2 = float(l.constant_term)
        x = (D*k1-B*k2)/(A*D-B*C)
        y = (A*k2-C*k1)/(A*D-B*C)
        print("The two lines intersects at x =", x, "and y =", y)
        return


    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n.coordinates)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


'''
l1 = Line(Vector([4.046,2.836]),1.21)
l2 = Line(Vector([10.115,7.09]),3.025)
l3 = Line(Vector([7.204,3.182]),8.68)
l4 = Line(Vector([8.172,4.114]),9.883)
l5 = Line(Vector([1.182,5.562]),6.744)
l6 = Line(Vector([1.773,8.343]),9.525)

l1.intersection(l2)
l3.intersection(l4)
l5.intersection(l6)

'''
