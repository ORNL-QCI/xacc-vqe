import pyxaccvqe as vqe
import pyxacc as xacc
from pyxaccvqe import PauliOperator
import unittest

class EigenDiagonalizeTest(unittest.TestCase):
        
    def testExecuteDiagonalize(self):
        src = """&FCI NORB=  2,NELEC=  2,MS2= 0,
            ORBSYM=1,5,
            ISYM=0,
            /                         i   a   j   b
            0.6744931033260081E+00   1   1   1   1
            0.6634720448605567E+00   2   2   1   1
            0.6973979494693358E+00   2   2   2   2
            0.1812875358123322E+00   2   1   2   1
            -0.1252477303982147E+01   1   1   0   0
            -0.4759344611440753E+00   2   2   0   0
            0.7137758743754461E+00   0   0   0   0
            """
        op = vqe.compile(src)
        self.assertEqual(vqe.execute(op).results[0][1], -1.137270422069755)

if __name__ == '__main__':
    unittest.main()
