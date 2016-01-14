/*
 * Bayesian Composite Models
 * 
 * (c) 2015 Collaborative Drug Discovery, Inc.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package com.cdd.bayes;

import com.cdd.bayes.ui.*;
import com.cdd.bayes.util.*;

import java.util.*;
import java.io.*;

/*
	Command line entrypoint: either 
*/

public class Main
{
	public static void main(String[] argv)
	{
		if (argv.length > 0 && (argv[0].equals("-h") || argv[0].equals("--help")))
		{
			printHelp();
			return;
		}
		
		Session session = new Session();
		
		final int len = argv.length;
		boolean openWindow = false;
		for (int n = 0; n < len; n++)
		{
			try
			{
    			if (argv[n].equals("-t"))
    			{
    				n++;
    				for (; n < len; n++) 
    				{
    					if (argv[n].startsWith("-")) {n--; break;}
    					session.addFile(parseDataFile(argv[n], Session.FILE_TRAINING));
    				}
    			}
    			else if (argv[n].equals("-s"))
    			{
    				n++;
    				for (; n < len; n++)
    				{
    					if (argv[n].startsWith("-")) {n--; break;}
	    				session.addFile(parseDataFile(argv[n], Session.FILE_TESTING));
    				}
    			}
    			else if (argv[n].equals("-p"))
    			{
    				n++;
    				for (; n < len; n++) 
    				{
    					if (argv[n].startsWith("-")) {n--; break;}
	    				session.addFile(parseDataFile(argv[n], Session.FILE_PREDICTION));
    				}
    			}
    			else if (argv[n].equals("-o"))
    			{
    				n++;
    				for (; n < len; n++)
    				{
    					if (argv[n].startsWith("-")) {n--; break;}
	    				session.addFile(parseDataFile(argv[n], Session.FILE_OUTPUT));
    				}
    			}
    			else if (argv[n].equals("-f") && n + 1 < len)
    			{
    				n++;
    				session.setFraction(Float.valueOf(argv[n]));
    			}
    			else if (argv[n].equals("-w")) openWindow = true;
    			else throw new IOException("Unexpected parameter.");
    		}
			catch (Exception ex)
			{
				Util.writeln("Unable to parse at parameter '" + argv[n] + "': " + ex.getMessage());
				return;
			}
		}
		
		if (openWindow || session.numFiles() == 0)
		{
			MainApplication.templateSession = session;
			new MainApplication().exec(new String[0]);
		}
		else
		{
			Util.writeln("Composite Bayesian Modelling: Session Parameters");
			for (int n = 0; n < session.numFiles(); n++)
			{
				Session.DataFile df = session.getFile(n);
				
				String strType = df.type == Session.FILE_TRAINING ? "Training" :
							     df.type == Session.FILE_TESTING ? "Testing" : 
							     df.type == Session.FILE_PREDICTION ? "Prediction" :
							     df.type == Session.FILE_OUTPUT ? "Output" : "?";
				String strField = df.field == null ? "(unknown)" : df.field;
				Util.writeln("    " + strType + " [" + df.filename + "] Field:[" + strField + "]");
			}
			Util.writeln("Fraction of training partitioned to testing set: " + session.getFraction());
		}
	}
	
	private static Session.DataFile parseDataFile(String fn, int type)
	{
		Session.DataFile df = new Session.DataFile(fn, type, null);
		int colon = fn.indexOf(':');
		if (colon >= 0) {df.filename = fn.substring(0, colon); df.field = fn.substring(colon + 1);}
		return df;
	}
	
	private static void printHelp()
	{
		Util.writeln("Composite Bayesian Modelling Tools\n");
		
		Util.writeln("Command line syntax:");
		Util.writeln("    -t <training files...>    input files for training set");
		Util.writeln("    -s <testing files...>     input files for testing set");
		Util.writeln("    -p <predicting files...>  input files for predictions");
		Util.writeln("    -o <output file>          output file to write predictions to");
		Util.writeln("    -f <fraction>             fraction (0..1) of training -> testing");
		Util.writeln("    -w                        open a window: interactive mode");

		Util.writeln("\nFiles can optionally be specified as <filename>:<fieldname>.");
		Util.writeln("Example syntax:");
		Util.writeln("    -t experimental.sdf:IC50 -f 0.1 -p newmolecules.sdf -o predicted.sdf");
	}
}
