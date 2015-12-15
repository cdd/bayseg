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

import com.cdd.bayes.util.*;

import java.util.*;

/*
 * Holds meta-data about a session, which describes the process of loading some number of input files for
 * training/testing/prediction, and the state involved with model generation.
 */
public class Session
{
	public static final int FILE_TRAINING = 1; // the {molecule:value} pairs are designated for training the model
	public static final int FILE_TESTING = 2; // the {molecule:value} pairs are for evaluating the model
	public static final int FILE_PREDICTION = 3; // molecules only: predictions are made of their possible value
	public static final int FILE_OUTPUT = 4; // output file for emitting predictions

	public static final class DataFile
	{
		public String filename;
		public int type; // one of FILE_*
		public String field;
		
		public DataFile(String filename, int type, String field)
		{
			this.filename = filename;
			this.type = type;
			this.field = field;
		}
	}
	private List<DataFile> files = new ArrayList<>();
	
	private float fraction = 0;
	
	// ------------ public methods ------------
	
	public Session()
	{
	}
	
	// data files: these make up the input/output
	public int numFiles() {return files.size();}
	public DataFile getFile(int idx) {return files.get(idx);}
	public void setFile(int idx, DataFile file) {files.set(idx, file);}
	public void addFile(DataFile file) {files.add(file);}
	public void deleteFile(int idx) {files.remove(idx);}
	
	// fraction of training set to push into the testing set
	public float getFraction() {return fraction;}
	public void setFraction(float fraction) {this.fraction = fraction;}
}


