import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class ViterbiTest {
	
    // test probOfSeq with example given in notes
	    @Test
	    public void probOfSeqExample() {
	        double[] pi = {.45, .55};
	        double[][] A = {{.7, .3},{.25,.75}};
	        double[][] B = {{.8, .2}, {.3, .7}};
		
	        Viterbi test = new Viterbi(pi, A, B);
		
	        int[] omega = {0,0,0,1};
			
	        double prob = test.probOfSeq(omega, test.bestSeq(omega, false), false);
	        assertTrue(Math.abs((prob - .0237)) < .0001);
	    }
	
	    // test that the maxScore of the sequence given in question 2.6 produces the given answer
	    @Test
	    public void question2_6() {
	        double[] pi = {.45, .55};
	        double[][] A = {{.7, .3},{.25,.75}};
	        double[][] B = {{.8, .2}, {.3, .7}};
		
	        Viterbi test = new Viterbi(pi, A, B);

	        int[] omega = new int[2000];
	        for (int i = 0; i < 500; i++) omega[i] = 0;
	        for (int i = 500; i < 1000; i++) omega[i] = 1;
	        for (int i = 1000; i < 1500; i++) omega[i] = 0;
	        for (int i = 1500; i < 2000; i++) omega[i] = 1;
		
	        double maxScore = test.maxScore(omega, true);
	        assertTrue(Math.abs((maxScore + 1227.5)) < .5);
	    }
	
	    // test that the maxScore of the sequence given in question 2.6 produces the given answer
		    @Test
		    public void question2_7() {
		        double[] pi = {.45, .55};
		        double[][] A = {{.7, .3},{.25,.75}};
		        double[][] B = {{.8, .2}, {.3, .7}};
			
		        Viterbi test = new Viterbi(pi, A, B);

		        int[] omega = new int[2004];
		        for (int i = 0; i < 500; i++) omega[i] = 0;
		        for (int i = 500; i < 1000; i++) omega[i] = 1;
		        for (int i = 1000; i < 1500; i++) omega[i] = 0;
		        for (int i = 1500; i < 2000; i++) omega[i] = 1;
		        omega[2000] = 0;
		        omega[2001] = 0;
		        omega[2002] = 0;
		        omega[2003] = 1;
			
		        double maxScore = test.maxScore(omega, true);
		        assertTrue(Math.abs((maxScore + 1231.8)) < .5);
		    }
	
		    // For the DNA example in question 3, test that maxScore is extremely close
		    // to the approximation given in the hw
		    @Test
		    public void dna_test() {
		        double[] pi_dna = {.5, .5};
		        double[][] A_dna = {{.5, .5},{.4,.6}};
		        double[][] B_dna = {{.2, .3, .3, .2}, {.3, .2, .2, .3}};
		
		        Viterbi dna = new Viterbi(pi_dna, A_dna, B_dna);
		        int[] omega_dna = {2,2,1,0,1,3,2,0,0};

		        List<Integer> list = new ArrayList<Integer>();
		        for (int x = 0; x < 3; x++) {
		            list.add(0);
		        }
		        for (int x = 0; x < 6; x++) {
		            list.add(1);
		        }
		
		        assertEquals(dna.bestSeq(omega_dna, false), list);

		        double maxScore = dna.maxScore(omega_dna, false);
		        // deal with approximation
		        assertTrue(Math.abs((maxScore - .00000004251)) < .00000000001);
		    }
	
		    // test that the maxScore of the sequence given in question 2.6 produces the given answer
		    @Test
		    public void underflow() {
		        double[] pi = {.45, .55};
		        double[][] A = {{.7, .3},{.25,.75}};
		        double[][] B = {{.8, .2}, {.3, .7}};
		
		        Viterbi test = new Viterbi(pi, A, B);

		        int[] omega = new int[1200];
		        for (int i = 0; i < 300; i++) omega[i] = 0;
		        for (int i = 300; i < 600; i++) omega[i] = 1;
		        for (int i = 600; i < 900; i++) omega[i] = 0;
		        for (int i = 900; i < 1200; i++) omega[i] = 1;
			
		        double maxScoreU = test.maxScore(omega, true);
		        double maxScore = test.maxScore(omega, false);
		        assertTrue(Math.exp(maxScoreU) == maxScore);
		    }
	
		    // Test that input edge cases throw the necessary exceptions
		    @Test
		    public void inputs() {
		        double[] pi = {.45, .55};
		        double[][] A = {{.7, .3},{.25,.75}};
		        double[][] B = {{.8, .2}, {.3, .7}};
		
		        Viterbi test = new Viterbi(pi, A, B);

		        int[] omega = null;
		
		        try {
		            test.checkOmega(omega);
		            fail("checkOmega with null input should throw exception!");
		        } catch (IllegalArgumentException e) {
			
		        }	
		        omega = new int[0];
		        try {
		            test.checkOmega(omega);
		            fail("checkOmega with input of 0 length should throw exception!");
		        } catch (IllegalArgumentException e) {
			
		        }
		
		        int k = -1;
		        try {
		            test.maxScoreSet(omega, k, false);
		            fail("if k is negative an exception should be thrown!");
		        } catch (IllegalArgumentException e) {
			
		        }
		        k = 10;
		        try {
		            test.maxScoreSet(omega, k, false);
		            fail("if k is greater than the number of distinct probabilities in the lattice an exception should be"
		                    + "thrown!");
		        } catch (IllegalArgumentException e) {
			
		        }
		
		        List<Integer> seq = null;
		        try {
		            test.probOfSeq(omega, seq, false);
		            fail("input of null sequence should throw exception");
		        } catch (IllegalArgumentException e) {
			
		        }
		        seq = new ArrayList<Integer>();
		        try {
		            test.probOfSeq(omega, seq, false);
		            fail("input of sequence of length 0 should throw exception");
		        } catch (IllegalArgumentException e) {
			
		        }
		
		    }
	
		    // tests some large inputs on functionality and performance
		    @Test
		    public void large() {
		        double[] pi = new double[500];
		        double[][] A = new double[500][500];
		        double[][] B = new double[500][500];
		        int[] omega = new int[500];
		        List<Integer> seq = new ArrayList<Integer>();
		
		        for (int i = 0; i < 500; i++) {
		            pi[i] = .002;
		            omega[i] = 0;
		            seq.add(0);
		            for (int j = 0; j < 500; j++) {
		                A[i][j] = .002;
		                B[i][j] = .002;
		            }
		        }
		        Viterbi test = new Viterbi(pi, A, B);
		
		        assertTrue(Math.abs(test.probOfSeq(omega, seq, true) + 6214.60809) < 1E-3);
		        List<Integer> bestSeq = test.bestSeq(omega, true);
		        List<Integer> worstSeq = test.worstSeq(omega, true);
		        List<List<Integer>> bestSeqSet = test.bestSeqSet(omega, true);
		        for (int i = 0; i < 500; i++) {
		            // out should be all 1's as each sequence has the same probability, so bestSeq should return
		            // the leftmost indexed state at each point
		            assertTrue(bestSeq.get(i) == 0);
		            assertTrue(worstSeq.get(i) == 0);
		            assertTrue(bestSeqSet.get(0).get(i) == 0);
		            assertTrue(bestSeqSet.get(1).get(i) == 499);
		        }
		
		        double maxScore = test.maxScore(omega, true);
		        double[] maxScoreSet = test.maxScoreSet(omega, 1, true);
		
		        assertTrue(Math.abs(maxScore + 6214.60809) < 1E-3);
		        assertTrue(Math.abs(maxScoreSet[0] + 6214.60809) < 1E-3);
		
		        try {
		            @SuppressWarnings("unused")
		            double[] maxScoreSetInvalid = test.maxScoreSet(omega, 10, true);
		            fail("maxScoreSetInvalid has < k distinct probabilities!");
		        } catch (IllegalArgumentException e) {
			
		        }
		    }

}
