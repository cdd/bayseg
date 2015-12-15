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

package com.cdd.bayes.ui;

import com.cdd.bayes.*;
import com.cdd.bayes.util.*;

import java.io.*;
import java.util.*;

import javafx.application.*;
import javafx.stage.*;
import javafx.scene.image.*;

/*
	BioAssay Ontology Tools: entrypoint with command line parameters.
*/

public class MainApplication extends Application
{
	public static Session templateSession = null; // typically constructed from the command line
	//public static Image icon = null;

	// ------------ public methods ------------	

	public MainApplication()
	{
		/*try
		{
    		InputStream istr = Util.openResource(this, "/images/MainIcon.png");
    		icon = new Image(istr);
    		istr.close();
		}
		catch (Exception ex) {ex.printStackTrace();}*/
	}
	
	public void exec(String[] args)
	{
		Application.launch(MainApplication.class, args);
	}
	
	public void start(Stage primaryStage)
	{
		CompositeWindow cw = new CompositeWindow(primaryStage, templateSession);
		final Stage stage = primaryStage;
        Platform.runLater(() -> stage.show());
	}

	// ------------ private methods ------------	

}
