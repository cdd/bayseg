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

import java.util.*;
import java.io.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.DefaultChemObjectBuilder;

/*
 * Takes a session instance as its parameter, and uses it to carry out the various steps involved with creating,
 * evaluating and applying a composite Bayesian operation. Some of the state is stored by modifying the session
 * object, while other parts are stored internally.
 */
public class ExecuteSession
{
	private Session session;
	private List<CompositeModel.Entry> training = new ArrayList<>(), testing = new ArrayList<>(), prediction = new ArrayList<>();
	private CompositeModel model = null;

	// ------------ public methods ------------
	
	public ExecuteSession(Session session)
	{
		this.session = session;
	}
	
	public Session getSession() {return session;}
	public List<CompositeModel.Entry> getTraining() {return training;}
	public List<CompositeModel.Entry> getTesting() {return testing;}
	public List<CompositeModel.Entry> getPrediction() {return prediction;}
	public CompositeModel getModel() {return model;}
	
	// loads the file indicated at the given index; clears out the previous batch of molecules; may fail gracefully (nop) or
	// complain with an exception
	public void loadFile(int idx) throws IOException
	{
		Session.DataFile df = session.getFile(idx);
		df.molecules.clear();
		if (df.filename == null || df.filename.length() == 0) return;
		File f = new File(df.filename);
		if (!f.exists()) throw new IOException("File not found: " + df.filename);
		if (!f.canRead()) throw new IOException("Access denied: " + df.filename);
		
        SDFixerHack hack = new SDFixerHack(new BufferedReader(new FileReader(f)));
        IteratingSDFReader rdr = new IteratingSDFReader(hack, DefaultChemObjectBuilder.getInstance());
        while (rdr.hasNext()) df.molecules.add(rdr.next());
        rdr.close();
	}
	
	// spool the loaded datafiles into the respective three partitions
	public void partitionMolecules()
	{
		training.clear();
		testing.clear();
		prediction.clear();
		
		for (Session.DataFile df : session.fileIter())
		{
			if (df.type == Session.FILE_OUTPUT) continue;
			if (df.type != Session.FILE_PREDICTION && df.field.length() == 0) continue;
			
			for (IAtomContainer mol : df.molecules)
			{
				CompositeModel.Entry entry = parseEntry(mol, df.type, df.field);
				if (entry == null) continue;
				if (df.type == Session.FILE_TRAINING) training.add(entry);
				else if (df.type == Session.FILE_TESTING) testing.add(entry);
				else if (df.type == Session.FILE_PREDICTION) prediction.add(entry);
			}
		}
		
		// if necessary, push some entries from training to testing
		int toMove = (int)Math.round(session.getFraction() * training.size());
		if (toMove > 0)
		{
			Random rnd = new Random(1); // predictable random
			while (toMove > 0 && training.size() > 10)
			{
				CompositeModel.Entry entry = training.remove(rnd.nextInt(training.size()));
				testing.add(entry);
				toMove--;
			}
		}
	}
	
	// stuff all the training set entries into the model and build it
	public void buildModel() throws CDKException
	{
		model = new CompositeModel();
		for (CompositeModel.Entry e : training) model.addEntry(e);
		model.determineSegments();
		model.calculate();
	}
	
	// ------------ private methods ------------

	// given a molecule that may or may not have an accompanying field datum, returns an entry: or null if not able to get enough
	// information out of it
	private CompositeModel.Entry parseEntry(IAtomContainer mol, int type, String field)
	{
		CompositeModel.Entry entry = new CompositeModel.Entry();
		entry.mol = mol;
		if (type == Session.FILE_PREDICTION) return entry;
		
		Object obj = mol.getProperties().get(field);
		if (obj == null || !(obj instanceof String)) return null;
		String str = (String)obj;

        while (str.startsWith(">") || str.startsWith("<") || str.startsWith(" ")) str = str.substring(1);
        try {entry.val = Double.parseDouble(str);}
        catch (NumberFormatException ex) {return null;}
        return entry;
	}

	/*
	 * An input stream intermediary which worksaround an unfortunate shortcoming in the CDK SD reader. For fields like:
	 * 
	 *     > <activity>
	 *     > 10
	 *     
	 * this is interpreted as two field declarations, so the value is skipped. This intermediary will convert such instances
	 * into:
	 * 
	 *     > <activity>
	 *     >10
	 *     
	 * which will be parsed as was originally intended.
	 */
	public final static class SDFixerHack extends InputStream
	{
		private BufferedReader in;
		private int pos = 0;
		private byte[] buff = null;

		public SDFixerHack(BufferedReader in)
		{
			this.in = in;
		}

		public int read() throws IOException
		{
			if (buff == null || pos >= buff.length) grabNextLine();
			if (buff == null) return -1;
			return buff[pos++];
		}

		public void close() throws IOException
		{
			in.close();
		}

		// perform the translation, if necessary
		private void grabNextLine() throws IOException
		{
			String line = in.readLine();
			if (line == null)
			{
				buff = null;
				return;
			}
			if (line.length() >= 3 && line.charAt(0) == '>' && line.charAt(1) == ' ' && line.charAt(2) != '<')
			{
				line = ">" + line.substring(2);
			}
			buff = (line + "\n").getBytes();
			pos = 0;
		}
	}	
}
