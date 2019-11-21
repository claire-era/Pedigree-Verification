package pedigreeVerification;

import java.util.ArrayList;

public class Main{
	public static void main(String[] args) {
		PreprocessHMP n = new PreprocessHMP(args);
		if(n.getVal()) {
			boolean eval = false;
			Step1 method1 = new Step1(args);
			int[] ascore = method1.getAutomatedScores();
			ArrayList<Integer> mscore = n.getManualScores();

			float count_1 = 0;
			
			for(int i = 0; i < mscore.size(); i++) {
				System.out.print(ascore[i] + ":" + mscore.get(i));
				System.out.print("\n");
				if(ascore[i] == 1 && mscore.get(i) == 1) {
					count_1++;
				}
			}
			
			//count yung mga number 1
			int count_ng_1 = 0;
			for(int j = 0; j < mscore.size(); j++) {
				if(mscore.get(j) == 1) count_ng_1++; 
			}
			Double rate = (double) (count_1 / count_ng_1);
			System.out.println((int) (rate * 100 ) + "% match");
		}else {
			System.out.println("may error sa data");
		}
	}
}