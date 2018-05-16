import pyxaccvqe as vqe
import pyxacc as xacc
from pyxaccvqe import PauliOperator
import unittest

class CompileTest(unittest.TestCase):
    def testFermionCompiler(self):
        
        src = """__qpu__ kernel() {
   0.7137758743754461
   -1.252477303982147 0 1 0 0
   0.337246551663004 0 1 1 1 1 0 0 0
   0.0906437679061661 0 1 1 1 3 0 2 0
   0.0906437679061661 0 1 2 1 0 0 2 0
   0.3317360224302783 0 1 2 1 2 0 0 0
   0.0906437679061661 0 1 3 1 1 0 2 0
   0.3317360224302783 0 1 3 1 3 0 0 0
   0.337246551663004 1 1 0 1 0 0 1 0
   0.0906437679061661 1 1 0 1 2 0 3 0
   -1.252477303982147 1 1 1 0
   0.0906437679061661 1 1 2 1 0 0 3 0
   0.3317360224302783 1 1 2 1 2 0 1 0
   0.0906437679061661 1 1 3 1 1 0 3 0
   0.3317360224302783 1 1 3 1 3 0 1 0
   0.3317360224302783 2 1 0 1 0 0 2 0
   0.0906437679061661 2 1 0 1 2 0 0 0
   0.3317360224302783 2 1 1 1 1 0 2 0
   0.0906437679061661 2 1 1 1 3 0 0 0
   -0.4759344611440753 2 1 2 0
   0.0906437679061661 2 1 3 1 1 0 0 0
   0.3486989747346679 2 1 3 1 3 0 2 0
   0.3317360224302783 3 1 0 1 0 0 3 0
   0.0906437679061661 3 1 0 1 2 0 1 0
   0.3317360224302783 3 1 1 1 1 0 3 0
   0.0906437679061661 3 1 1 1 3 0 1 0
   0.0906437679061661 3 1 2 1 0 0 1 0
   0.3486989747346679 3 1 2 1 2 0 3 0
   -0.4759344611440753 3 1 3 0
}"""

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
