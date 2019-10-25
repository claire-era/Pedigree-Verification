package rtests;

public class Main {
	public static void main(String[] args) {	
		rtest a = new rtest();
		Double pval = a.compute();
		Double pval2 = a.compute();
		System.out.println(pval + " " + pval2);
	}
}
