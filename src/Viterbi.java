import static org.junit.Assert.assertTrue;
import java.util.ArrayList;
import java.util.List;

public class Viterbi {
    
    double[] pi;
	    double[][] A;
	    double[][] B;
	    // used to avoid printing unnecessarily when a function calls another one
	    boolean print = true;
	
	    // creates a new Viterbi-type gsm  
	    public Viterbi (double[] _pi, double[][] _A, double[][] _B) {
	        pi = _pi;
	        A = _A;
	        B = _B;
		
	        // Check that pi represents a proper vector of probabilities and that dimensions are acceptable
	        if (pi.length != A.length) {
             throw new java.lang.IllegalArgumentException("Input matrix pi does not match the dimensions of matrix A!");
	        }
	        double pi_sum = 0;
	        for (int i = 0; i < pi.length; i++) {
	            pi_sum += pi[i];
	        }
	        if (Math.abs(pi_sum - 1) > 1E-15) {
	            System.out.println("qwuodb ---" + pi_sum);
	            throw new java.lang.IllegalArgumentException(
	                    "Input matrix pi does not represent a proper vector of probabilities!");
	        }
	
	        // Check that A represents a proper matrix of probabilities and that dimensions are acceptable
	        for (int i = 0; i < A.length; i++) {
	            if (A.length != A[i].length) {
	                throw new java.lang.IllegalArgumentException("Input matrix A has malformed dimensions!");
	            }
			
	            double sum = 0;
	            for (int j = 0; j < A.length; j++) {
	                sum += A[i][j];
	            }
	            if (Math.abs(sum - 1) > 1E-15) {
	                throw new java.lang.IllegalArgumentException(
	                        "Input matrix A does not represent a proper matrix of probabilities!");
	            }
	        }
		
	        // Check that B represents a proper matrix of probabilities and that dimensions are acceptable
	        for (int i = 0; i < B.length; i++) {
	            if (B[i].length != B[0].length) {
	                throw new java.lang.IllegalArgumentException("Input matrix B has malformed dimensions!");
	            }
			
	            double sum = 0;
	            for (int j = 0; j < B[0].length; j++) {
	                sum += B[i][j];
	            }
	            if (Math.abs(sum - 1) > 1E-15) {
	                throw new java.lang.IllegalArgumentException(
	                        "Input matrix B does not represent a proper matrix of probabilities!");
	            }
	        }
		
	    }

	
	    // computes and prints the probability of a particular sequence of occurring
	    public double probOfSeq(int[] omega, List<Integer> seq, boolean underflow) {
	        checkOmega(omega);
	        if (seq == null) {
	            throw new java.lang.IllegalArgumentException("seq is null!");
	        } else if (seq.size() == 0) {
	            throw new java.lang.IllegalArgumentException("seq size is zero!");
	        }
		
	        double prob = pi[seq.get(0)] * B[seq.get(0)][omega[0]];
	        if (underflow) {
	            prob = log(pi[seq.get(0)]) + log(B[seq.get(0)][omega[0]]);
	        }
	        if (omega.length > 1) {
	            for (int t = 1; t < seq.size(); t++) {
	                if (!underflow) {
	                    prob *= A[seq.get(t-1)][seq.get(t)]*B[seq.get(t)][omega[t]]; 
	                } else {
	                    prob += log(A[seq.get(t-1)][seq.get(t)]) + log(B[seq.get(t)][omega[t]]);
	                }
	            }
	        }

	        System.out.println(prob);
	        return prob;
	    }
	
	
	    // computes and prints the sequence with the highest probability of occurring
	    public List<Integer> bestSeq(int[] omega, boolean underflow) {
	        checkOmega(omega);

	        double[][] scores = new double[A.length][omega.length];
		
	        // initial scores based on start state probabilities
	        for (int j = 0; j < pi.length; j++) {
	            if (!underflow) {
	                scores[j][0] = pi[j]*B[j][omega[0]];
	            } else {
	                scores[j][0] = log(pi[j]) + log(B[j][omega[0]]);
	            }
	        }
		
	        // if omega is of length one, we can output the most likely state that produced the single output
	        if (omega.length == 1) {
	            List<Integer> out = new ArrayList<Integer>();
	            int maxI = 0;
	            for (int j = 0; j < pi.length; j++) {
	                if (scores[j][0] > scores[maxI][0]) {
	                    maxI = j;
	                }
	            }
	            out.add(maxI);
	            System.out.println(maxI);
	            return out;
	        }
		
	        // Stores all pred[i][t], which is equal to the state preceding i in the best sequence of length t
	        int[][] pred = new int[A.length][omega.length];
	
	        for (int t = 1; t < omega.length; t++) {
			
	            for (int j = 0; j < A.length; j++) {
	                double maxTScore = 0;
	                if (underflow) {
	                    maxTScore = Double.NEGATIVE_INFINITY;
	                }
	                int maxI = 0;
				
	                for (int i = 0; i < A.length; i++) {
	                    double temp = scores[i][t-1]*A[i][j]*B[j][omega[t]];
	                    if (underflow) {
	                        temp = scores[i][t-1] + log(A[i][j]) + log(B[j][omega[t]]);
	                    }
	                    if (temp > maxTScore) {
	                        maxTScore = temp;
	                        maxI = i;
	                    }
	                }
	                scores[j][t] = maxTScore;
	                pred[j][t] = maxI;
	            }
	        }	

	        List<Integer> out = new ArrayList<Integer>();

	        int lastState = 0;
	        double max = 0;
	        if (underflow) {
	            max = Double.NEGATIVE_INFINITY;
	        }
	        for (int y = 0; y < A.length; y++) {
	            if (scores[y][omega.length-1] > max) {
	                max = scores[y][omega.length-1];
	                lastState = y;
	            }
	        }
	        out.add(lastState);
		
	        for (int x = omega.length-1; x > 0; x--) {
	            lastState = pred[lastState][x];
	            out.add(0, lastState);
	        }
		
	        if (print) {
	            printList(out);
	        }
	        return out;
	    }

	
	    // computes and prints the sequence with the lowest probability of occurring
	    public List<Integer> worstSeq(int[] omega, boolean underflow) {
	        checkOmega(omega);

	        double[][] scores = new double[A.length][omega.length];
		
	        // initial scores based on start state probabilities
	        for (int j = 0; j < pi.length; j++) {
	            if (!underflow) {
	                scores[j][0] = pi[j]*B[j][omega[0]];
	            } else {
	                scores[j][0] = log(pi[j]) + log(B[j][omega[0]]);
	            }
	        }
		
	        // if omega is of length one, we can output the most likely state that produced the single output
	        if (omega.length == 1) {
	            List<Integer> out = new ArrayList<Integer>();
	            int minI = 0;
	            for (int j = 0; j < pi.length; j++) {
	                if (scores[j][0] < scores[minI][0]) {
	                    minI = j;
	                }
	            }
	            out.add(minI);
	            System.out.println(minI);
	            return out;
	        }
		
	        // Stores all pred[i][t], which is equal to the state preceding i in the best sequence of length t
	        int[][] pred = new int[A.length][omega.length];
	
	        for (int t = 1; t < omega.length; t++) {
			
	            for (int j = 0; j < A.length; j++) {
	                double minTScore = 1;
	                int minI = 0;
				
	                for (int i = 0; i < A.length; i++) {
	                    double temp = scores[i][t-1]*A[i][j]*B[j][omega[t]];
	                    if (underflow) {
	                        temp = scores[i][t-1] + log(A[i][j]) + log(B[j][omega[t]]);
	                    }
	                    if (temp < minTScore) {
	                        minTScore = temp;
	                        minI = i;
	                    }
	                }
	                scores[j][t] = minTScore;
	                pred[j][t] = minI;
	            }
	        }	

	        List<Integer> out = new ArrayList<Integer>();

	        int lastState = 0;
	        double min = 1;
	        for (int y = 0; y < A.length; y++) {
	            if (scores[y][omega.length-1] < min) {
	                min = scores[y][omega.length-1];
	                lastState = y;
	            }
	        }
	        out.add(lastState);
		
	        for (int x = omega.length-1; x > 0; x--) {
	            lastState = pred[lastState][x];
	            out.add(0, lastState);
	        }
		
	        printList(out);
	        return out;
		
	    }
	
	
	    // for t = 1, …, T and j = 1, …, n, compute predl(j,t) to be the lowest index l
	    // such that tscore(l) (at time t) is maximum, and predr(j, t) to be the largest
	    // index r such that tscore(r) (at time t) is maximum. If l = r, then predl(j,t) =
	    // predr(j,t). Compute and print the sequence with the highest probability
	    // corresponding to the indices in predl, and the sequence with the highest
	    // probability corresponding to the indices in predr. (12​ ​points) 
	    public List<List<Integer>> bestSeqSet(int[] omega, boolean underflow) {
	        checkOmega(omega);

	        double[][] scores = new double[A.length][omega.length];
		
	        // initial scores based on start state probabilities
	        for (int j = 0; j < pi.length; j++) {
	            if (!underflow) {
	                scores[j][0] = pi[j]*B[j][omega[0]];
	            } else {
	                scores[j][0] = log(pi[j]) + log(B[j][omega[0]]);
	            }
	        }
		
	        // if omega is of length one, we can output the most likely state that produced the single output
	        if (omega.length == 1) {
	            List<List<Integer>> out = new ArrayList<List<Integer>>();
			
	            List<Integer> maxI = new ArrayList<Integer>();
	            for (int j = 0; j < pi.length; j++) {
	                if (scores[j][0] > scores[maxI.get(0)][0]) {
	                    maxI = new ArrayList<Integer>();
	                    maxI.add(j);
	                } else if (scores[j][0] == scores[maxI.get(0)][0]) {
	                    maxI.add(j);
	                }
	            }
	            out.add(new ArrayList<Integer>(maxI.get(0)));
	            out.add(new ArrayList<Integer>(maxI.get(maxI.size()-1)));
	            for (List<Integer> el : out) {
	                printList(el);
	            }
	            return out;
	        }
		
	        // Stores all predl[i][t], which is the lowest index l such that tscore(l) (at time t) is maximum
	        int[][] predl = new int[A.length][omega.length];
	        // Stores all predr[i][t], which is the largest index r such that tscore(r) (at time t) is maximum
	        int[][] predr = new int[A.length][omega.length];
	
	        for (int t = 1; t < omega.length; t++) {
			
	            for (int j = 0; j < A.length; j++) {
	                double maxTScore = 0;
	                if (underflow) {
	                    maxTScore = Double.NEGATIVE_INFINITY;
	                }
	                int maxIL = 0;
	                int maxIR = 0;
				
	                for (int i = 0; i < A.length; i++) {
	                    double temp = scores[i][t-1]*A[i][j]*B[j][omega[t]];
	                    if (underflow) {
	                        temp = scores[i][t-1] + log(A[i][j]) + log(B[j][omega[t]]);
	                    }
	                    if (temp > maxTScore) {
	                        maxTScore = temp;
	                        maxIL = i;
	                        maxIR = i;
	                    } else if (temp == maxTScore) {
	                        maxIR = i;
	                    }
	                }
	                scores[j][t] = maxTScore;
	                predl[j][t] = maxIL;
	                predr[j][t] = maxIR;
	            }
	        }	

	        // Output is a list of lists, the leftmost highest score path, an the rightmost one
	        List<List<Integer>> out = new ArrayList<List<Integer>>();
	        List<Integer> l = new ArrayList<Integer>();
	        List<Integer> r = new ArrayList<Integer>();
		
	        int lastStateL = 0;
	        int lastStateR = 0;
	        double max = 0;
	        if (underflow) {
	            max = Double.NEGATIVE_INFINITY;
	        }
	        for (int y = 0; y < A.length; y++) {
	            if (scores[y][omega.length-1] > max) {
	                max = scores[y][omega.length-1];
	                lastStateL = y;
	                lastStateR = y;
	            } else if (scores[y][omega.length-1] == max) {
	                lastStateR = y;
	            }
	        }
	        l.add(lastStateL);
	        r.add(lastStateR);
		
	        for (int x = omega.length-1; x > 0; x--) {
	            lastStateL = predl[lastStateL][x];
	            l.add(0, lastStateL);
			
	            lastStateR = predr[lastStateR][x];
	            r.add(0, lastStateR);
	        }
		
	        out.add(l);
	        out.add(r);

	        for (List<Integer> el : out) {
	            printList(el);
	        }
	        return out;
	    }
	
	
	    // computes and prints the highest probability found at time T. 
	    // (find probability of best sequence)
	    public double maxScore(int[] omega, boolean underflow) {
	        // print just stops bestSeq() from printing output
	        print = false;
	        List<Integer> seq = bestSeq(omega, underflow);
	        print = true;
		
	        return probOfSeq(omega, seq, underflow);
		
	    }
	
	
	    // computes and prints the k distinct highest probabilities found at time T
	    // (find the k probabilities of the k best sequences)
	    public double[] maxScoreSet(int[] omega, int k, boolean underflow) {
	        checkOmega(omega);
	        if (k > A.length) {
	            throw new IllegalArgumentException("k is larger than the number of distinct probabilities in the lattice!");
	        } else if (k < 0) {
	            throw new IllegalArgumentException("k is less than or equal to 0!");
	        }
		
	        double[][] scores = new double[A.length][omega.length];
		
	        // initial scores based on start state probabilities
	        for (int j = 0; j < pi.length; j++) {
	            if (!underflow) {
	                scores[j][0] = pi[j]*B[j][omega[0]];
	            } else {
	                scores[j][0] = log(pi[j]) + log(B[j][omega[0]]);
	            }
	        }
		
	        // if omega is of length one, we can output the most likely state that produced the single output
	        if (omega.length == 1) {
	            double[] out = new double[1];
	            int maxI = 0;
	            for (int j = 0; j < pi.length; j++) {
	                if (scores[j][0] > scores[maxI][0]) {
	                    maxI = j;
	                }
	            }
	            out[0] = scores[maxI][0];
	            System.out.println(out[0]);
	            return out;
	        }
		
	        // Stores all pred[i][t], which is equal to the state preceding i in the best sequence of length t
	        int[][] pred = new int[A.length][omega.length];
	
	        for (int t = 1; t < omega.length; t++) {
			
	            for (int j = 0; j < A.length; j++) {			
	                double maxTScore = 0;
	                if (underflow) {
	                    maxTScore = Double.NEGATIVE_INFINITY;
	                }
	                int maxI = 0;
				
	                for (int i = 0; i < A.length; i++) {
	                    double temp = scores[i][t-1]*A[i][j]*B[j][omega[t]];
	                    if (underflow) {
	                        temp = scores[i][t-1] + log(A[i][j]) + log(B[j][omega[t]]);
	                    }
	                    if (temp > maxTScore) {
	                        maxTScore = temp;
	                        maxI = i;
	                    }
	                }
	                scores[j][t] = maxTScore;
	                pred[j][t] = maxI;
	            }
	        }	
		
	        double[] out = new double[k];
	        for (int x = 0; x < k; x++) {
	            double max = 0;
	            if (underflow) {
	                max = Double.NEGATIVE_INFINITY;
	            }
	            for (int y = 0; y < A.length; y++) {
	                if (scores[y][omega.length-1] > max) {
	                    boolean alreadyAdded = false;
	                    for (int z = 0; z < x; z++) {
	                        if (out[z] == scores[y][omega.length-1]) {
	                            alreadyAdded = true;
	                        }
	                    }
	                    if (!alreadyAdded) {
	                        max = scores[y][omega.length-1];
	                    }
	                }
	            }
	            if ((max == 0) || (max == Double.NEGATIVE_INFINITY)) {
	                throw new IllegalArgumentException("k is invalid, there aren't that many distinct probabilities!");
	            }
	            out[x] = max;
	        }
		
	        for (int x = 0; x < out.length; x++) {
	            System.out.println(out[x]);
	        }
	        return out;		
	    }
	
	
	    // Checks if omega is null or of length 0, then throws the necessary IllegalArgumentException
	    public void checkOmega(int[] omega) {
	        if (omega == null) {
	            throw new java.lang.IllegalArgumentException("omega is null!");
	        } else if (omega.length == 0) {
	            throw new java.lang.IllegalArgumentException("omega has 0 length!");
	        }
	    }
	
	    // Outputs log(input). This just uses the Math.log function, but checks for input of 0 and returns
	    // Double.NEGATIVE_INFINITY in this special case.
	    public double log(double a) {
	        if (a == 0) {
	            return Double.NEGATIVE_INFINITY;
	        }
	        return Math.log(a);
	    }
	
	    public void printList(List<Integer> list) {
	        if (list.size() == 0) {
	            System.out.println();
	        } else {
	            System.out.print(list.get(0));
	            for (int i = 1; i < list.size(); i++) {
	                System.out.print("," + list.get(i));
	            }
	            System.out.println();
	        }		
	    }
	
}



