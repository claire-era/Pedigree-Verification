package rtests;
import java.io.*;
import java.awt.Frame;
import java.awt.FileDialog;

import java.util.Enumeration;

import org.rosuda.JRI.Rengine;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RList;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.RMainLoopCallbacks;

public class FisherTest {
    public static void main(String[] args) {
	// just making sure we have the right version of everything
    	System.out.println(System.getProperty("java.library.path"));
    	
	if (!Rengine.versionCheck()) {
	    System.err.println("** Version mismatch - Java files don't match library version.");
	    System.exit(1);
	}
        System.out.println("Creating Rengine (with arguments)");
		// 1) we pass the arguments from the command line
		// 2) we won't use the main loop at first, we'll start it later
		//    (that's the "false" as second argument)
		// 3) the callbacks are implemented by the TextConsole class above
		Rengine re=new Rengine(args, false, new TextConsole());
		// the engine creates R is a new thread, so we should wait until it's ready
        if (!re.waitForR()) {
            System.out.println("Cannot load R");
            return;
        }
        
		try {
//			REXP x;
//			REXP matrix;
			String sample_m = "matrix(c(3, 1, 1, 3), nrow = 2, dimnames = list(Guess = c(\"Milk\", \"Tea\"),Truth = c(\"Milk\", \"Tea\")))";
			re.eval("TeaTasting=" + sample_m);
			re.eval("result=fisher.test(TeaTasting, alternative = \"greater\")");
			Double pval = re.eval("result$p.value").asDouble();
			System.out.println("P-value is =" + pval);
		}finally {}
        
	// so far we used R as a computational slave without REPL
	// now we start the loop, so the user can use the console
//	System.out.println("Now the console is yours ... have fun");
//	re.startMainLoop();
    }
}
