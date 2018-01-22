import pyxaccvqe as vqe
import pyxacc as xacc
from pyxaccvqe import PauliOperator
import unittest

class CompileTest(unittest.TestCase):
    def testFCIDump(self):
        
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
        
        self.assertEqual(op.nTerms(), 15)
        
        expected = PauliOperator({2:'Z', 3:'Z'}, .174349)
        expected += PauliOperator({1:'Z', 3:'Z'}, .120546)
        expected += PauliOperator({3:'Z'}, -.222796)
        expected += PauliOperator({0:'Z', 3:'Z'}, .165868)
        expected += PauliOperator({0:'Z', 2:'Z'}, .120546)
        expected += PauliOperator({1:'Z', 2:'Z'}, .165868)
        expected += PauliOperator({0:'Z', 1:'Z'}, .168623)
        expected += PauliOperator({1:'Z'}, .171201)
        expected += PauliOperator({2:'Z'}, -.222796)
        expected += PauliOperator({0:'Z'}, .171201)
        expected += PauliOperator({}, -.0988349)
        expected += PauliOperator({0:'X', 1:'X', 2:'Y', 3:'Y'}, -0.0453219)
        expected += PauliOperator({0:'Y', 1:'X', 2:'X', 3:'Y'}, 0.0453219)
        expected += PauliOperator({0:'Y', 1:'Y', 2:'X', 3:'X'}, -0.0453219)
        expected += PauliOperator({0:'X', 1:'Y', 2:'Y', 3:'X'}, 0.0453219)
        
        self.assertTrue(op == expected)
        self.assertTrue(op.isClose(expected))
        
if __name__ == '__main__':
    unittest.main()
