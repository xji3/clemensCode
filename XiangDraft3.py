from protein import *
from data_test import *
from operator import attrgetter

class Student:
    def __init__(self, name, grade, age):
            self.name = name
            self.grade = grade
            self.age = age
    def __repr__(self):
            return repr((self.name, self.grade, self.age))


student_objects = [
        Student('john', 'C', 15),
        Student('jane', 'A', 12),
        Student('dave', 'B', 10),
        Student('Adam', 'A', 9)
]

print sorted(student_objects,key=attrgetter('age'))
print sorted(student_objects,key=attrgetter('grade'))
