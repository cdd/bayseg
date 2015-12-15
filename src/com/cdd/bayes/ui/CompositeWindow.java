/*
 * BioAssay Ontology Annotator Tools
 * 
 * (c) 2014-2015 Collaborative Drug Discovery Inc.
 */

package com.cdd.bayes.ui;

import com.cdd.bayes.*;
import com.cdd.bayes.util.*;

import java.io.*;
import java.net.*;
import java.util.*;

import javafx.event.*;
import javafx.geometry.*;
import javafx.stage.*;
import javafx.scene.*;
import javafx.scene.control.*;
import javafx.scene.input.*;
import javafx.scene.layout.*;
import javafx.scene.paint.*;
import javafx.application.*;
import javafx.beans.value.*;
import javafx.util.*;

import org.json.*;

/*
*/

public class CompositeWindow
{
	// ------------ private data ------------	

	private Session session;

    private Stage stage;
    private BorderPane root;
    
    //private MenuBar menuBar;
    //private Menu menuFile, menuEdit, menuValue, menuView;
    
	// ------------ public methods ------------	

	public CompositeWindow(Stage stage, Session session)
	{
		this.stage = stage;
		this.session = session;

        stage.setTitle("Composite Bayesian Modelling");

		/*menuBar = new MenuBar();
		menuBar.setUseSystemMenuBar(true);
		menuBar.getMenus().add(menuFile = new Menu("_File"));
		menuBar.getMenus().add(menuEdit = new Menu("_Edit"));
		menuBar.getMenus().add(menuValue = new Menu("_Value"));
		menuBar.getMenus().add(menuView = new Menu("Vie_w"));
		createMenuItems();*/

		root = new BorderPane();
		//root.setTop(menuBar);
		root.setCenter(new Label("fnord!"));

		Scene scene = new Scene(root, 700, 600, Color.WHITE);

		stage.setScene(scene);
		
		//new Thread(() -> backgroundLoadTemplates()).run();
 	}


	// ------------ private methods ------------	
}
