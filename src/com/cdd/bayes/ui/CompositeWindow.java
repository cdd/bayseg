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

import javafx.scene.control.*;
import javafx.event.*;
import javafx.geometry.*;
import javafx.stage.*;
import javafx.scene.*;
import javafx.scene.control.*;
import javafx.scene.input.*;
import javafx.scene.layout.*;
import javafx.scene.paint.*;
import javafx.scene.shape.*;
import javafx.application.*;
import javafx.beans.value.*;
import javafx.util.*;

import org.json.*;
import org.controlsfx.control.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.exception.CDKException;

/*
 * Interactive window display for setting up, executing and visualizing a composite Bayesian session.
*/

public class CompositeWindow
{
	// ------------ private data ------------	

	private Session session; // when "busy", this is owned by a background task; otherwise UI thread can use it

    private Stage stage;
    private VBox content = new VBox();
    
	final int PADDING = 4;

    //private MenuBar menuBar;
    //private Menu menuFile, menuEdit, menuValue, menuView;
    
    // widgets that get updated in sync with the session
    private static final class FileGroup
    {
    	VBox vbox;
    	Label labelFilename;
    	Button btnChoose, btnDelete;
    	SegmentedButton segType;
    	ToggleButton radioTrain, radioTest, radioPred, radioOutput;
    	ComboBox<String> comboField;
    }
    private List<FileGroup> fileGroups = new ArrayList<>();
    private TextField textFraction;
    private Button btnAddFile, btnLoad, btnBuild, btnPredict, btnSave;

    private Object mutex = new String("!");
    private boolean busy = false; // set to true if doing something behind the scenes; wrap with mutex
    private ExecuteSession exec = null; // as with session: busy = owned by background task
    
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

		ScrollPane scroller = new ScrollPane(content);
		scroller.setHbarPolicy(ScrollPane.ScrollBarPolicy.NEVER);
		scroller.setVbarPolicy(ScrollPane.ScrollBarPolicy.AS_NEEDED);
		scroller.setFitToWidth(true);
		content.setPadding(new Insets(PADDING));
		content.setPrefWidth(Double.MAX_VALUE);


		BorderPane root = new BorderPane();
		//root.setTop(menuBar);
		root.setCenter(scroller);

		Scene scene = new Scene(root, 700, 700, Color.WHITE);

		stage.setScene(scene);

		recreateContent();
 	}

	// ------------ private methods ------------	
	
	private void recreateContent()
	{
		content.getChildren().clear();
		fileGroups.clear();
		
		content.setSpacing(PADDING * 2);
		
		for (Session.DataFile df : session.fileIter()) 
		{
			final  int idx = fileGroups.size();
		
			FileGroup fg = new FileGroup();
			fg.vbox = new VBox();
			fg.vbox.setPrefWidth(Double.MAX_VALUE);
			
			fg.labelFilename = new Label(df.filename);
			fg.labelFilename.setStyle("-fx-text-fill: black; -fx-border-color: black; -fx-background-color: #F0F0F0; -fx-padding: 0.2em 0.2em 0.2em 0.2em;");
			fg.labelFilename.setTextOverrun(OverrunStyle.LEADING_ELLIPSIS);
			fg.labelFilename.setMinWidth(50);
			fg.labelFilename.setMaxWidth(Double.MAX_VALUE);
			fg.btnChoose = new Button("\u21b4");
			fg.btnChoose.setTooltip(new Tooltip("Select a different file"));
			fg.btnChoose.setOnAction((value) -> actionChooseFile(idx));
			fg.btnDelete = new Button("\u2716");
			fg.btnDelete.setTooltip(new Tooltip("Remove file"));
			fg.btnDelete.setOnAction((value) -> actionDeleteFile(idx));
			
			Lineup line = new Lineup(PADDING);

			RowLine row = new RowLine(PADDING);
			row.add(fg.labelFilename, 1);
			row.add(fg.btnChoose, 0);
			row.add(fg.btnDelete, 0);
			line.add(row, "File:", 1, 0);

			fg.radioTrain = new ToggleButton("Training");
			fg.radioTest = new ToggleButton("Testing");
			fg.radioPred = new ToggleButton("Prediction");
			fg.radioOutput = new ToggleButton("Output");
			fg.radioTrain.setOnAction((value) -> changeType(idx, Session.FILE_TRAINING));
			fg.radioTest.setOnAction((value) -> changeType(idx, Session.FILE_TESTING));
			fg.radioPred.setOnAction((value) -> changeType(idx, Session.FILE_PREDICTION));
			fg.radioOutput.setOnAction((value) -> changeType(idx, Session.FILE_OUTPUT));
			fg.radioTrain.setTooltip(new Tooltip("Molecules and activities used for training the models"));
			fg.radioTest.setTooltip(new Tooltip("Molecules and activities used for evaluating models"));
			fg.radioPred.setTooltip(new Tooltip("Molecules are used to make predictions"));
			fg.radioOutput.setTooltip(new Tooltip("Predictions will be written to this file"));
			fg.segType = new SegmentedButton(fg.radioTrain, fg.radioTest, fg.radioPred, fg.radioOutput);
			if (df.type == Session.FILE_TRAINING) fg.radioTrain.setSelected(true);
			else if (df.type == Session.FILE_TESTING) fg.radioTest.setSelected(true);
			else if (df.type == Session.FILE_PREDICTION) fg.radioPred.setSelected(true);
			else if (df.type == Session.FILE_OUTPUT) fg.radioOutput.setSelected(true);
			
			fg.comboField = new ComboBox<>();
			fg.comboField.setEditable(true);
			fg.comboField.setValue(df.field);
			fg.comboField.setMaxWidth(Double.MAX_VALUE);
			fg.comboField.setTooltip(new Tooltip("Specify the field used for activity data"));
			fg.comboField.valueProperty().addListener((observe, oldval, newval) -> changeField(idx, newval));
			fg.comboField.getEditor().textProperty().addListener((observe, oldval, newval) -> changeField(idx, newval));
			
			Lineup lineField = new Lineup(PADDING);
			lineField.add(fg.comboField, "Field:", 1, 0);
			
			row = new RowLine(PADDING);
			row.add(fg.segType, 0);
			int nmol = df.molecules.size();
			if (nmol > 0) 
			{
				Label labelRows = new Label("(" + nmol + " row" + (nmol == 1 ? "" : "s") + ")");
				labelRows.setStyle("-fx-font-style: italic; -fx-text-fill: #404040;");
				row.add(labelRows);

			}
			row.add(lineField, 1);
			
			line.add(row, "Type:", 1, 0);

			fg.vbox.getChildren().add(line);
			
			content.getChildren().add(fg.vbox);
			fileGroups.add(fg);

			// scan through the molecules and add all unique fields to the combobox			
			Set<String> already = new HashSet<>();
			for (IAtomContainer mol : df.molecules)
			{
				for (Object obj : mol.getProperties().keySet()) if (obj instanceof String)
				{
					String fldName = (String)obj;
					if (already.contains(fldName)) continue;
					fg.comboField.getItems().add(fldName);
					already.add(fldName);
				}
			}
			
			//content.getChildren().add(new Rectangle(0, 10)); // spacer
		}
		
		Lineup line = new Lineup(PADDING);
		textFraction = new TextField(String.valueOf(session.getFraction()));
		textFraction.textProperty().addListener((observe, oldval, newval) -> changeFraction(newval));		
		textFraction.setTooltip(new Tooltip("Fraction of training set set aside for testing."));
		line.add(textFraction, "Reserve Testing Fraction:", 0, 0);
		content.getChildren().add(line);

		content.getChildren().add(new Separator());

		btnAddFile = new Button("Add File");
		btnLoad = new Button("Load");
		btnBuild = new Button("Build");
		btnPredict = new Button("Predict");
		btnSave = new Button("Save");
		btnAddFile.setTooltip(new Tooltip("Select new file to use"));
		btnLoad.setTooltip(new Tooltip("Load all input datafiles"));
		btnBuild.setTooltip(new Tooltip("Build model from input files"));
		btnPredict.setTooltip(new Tooltip("Make predictions using model"));
		btnSave.setTooltip(new Tooltip("Save predictions to output file"));
		btnAddFile.setOnAction((value) -> actionAddFile());
		btnLoad.setOnAction((value) -> actionLoad());
		btnBuild.setOnAction((value) -> actionBuild());
		btnPredict.setOnAction((value) -> actionPredict());
		btnSave.setOnAction((value) -> actionSave());
		
		FlowPane flow = new FlowPane();
		flow.setAlignment(Pos.BASELINE_RIGHT);
		flow.setHgap(PADDING);
		flow.getChildren().addAll(btnAddFile, btnLoad, btnBuild, btnPredict, btnSave);
		content.getChildren().add(flow);
		
		if (exec != null) 
		{
			content.getChildren().add(new Separator());
			recreateExecution();
		}
		
		updateContent();
	}
	
	private void recreateExecution()
	{
		List<CompositeModel.Entry> training = exec.getTraining(), testing = exec.getTesting(), prediction = exec.getPrediction();
	
		if (training.size() > 0 || testing.size() > 0 || prediction.size() > 0)
		{
    		FlowPane flow = new FlowPane();
    		flow.setAlignment(Pos.BASELINE_LEFT);
    		flow.setHgap(PADDING);
		
			Label labelTrain = new Label(String.valueOf(exec.getTraining().size()));
			Label labelTest = new Label(String.valueOf(exec.getTesting().size()));
			Label labelPred = new Label(String.valueOf(exec.getPrediction().size()));
			for (Label label : new Label[]{labelTrain, labelTest, labelPred}) label.setStyle("-fx-text-fill: black; -fx-border-color: black; -fx-background-color: #F8F8F8; -fx-padding: 0.2em 0.2em 0.2em 0.2em;");
			flow.getChildren().addAll(new Label("Size of Training Set:"), labelTrain);
			flow.getChildren().addAll(new Label("Testing Set:"), labelTest);
			flow.getChildren().addAll(new Label("Prediction Set:"), labelPred);
			
			content.getChildren().add(flow);
		}
		
		CompositeModel model = exec.getModel();
		
		// show the data distribution for training & testing
		boolean invertDir = false;
		if (training.size() > 0 || testing.size() > 0)
		{
			double[] segments = model == null ? null : model.getSegments();
			
			List<List<CompositeModel.Entry>> datasets = new ArrayList<>();
			datasets.add(exec.getTraining());
			if (exec.getTesting().size() > 0) datasets.add(exec.getTesting());
		
			RenderDatasets render = new RenderDatasets(datasets, segments);
			render.draw();
			invertDir = render.usesLog();
			
			HBox hbox = new HBox();
			hbox.setAlignment(Pos.TOP_CENTER);
			hbox.getChildren().add(render.getCanvas());
			content.getChildren().add(hbox);
		}
		
		// show the cross validation matrices
		if (model != null)
		{
			HBox hbox = new HBox();
			hbox.setAlignment(Pos.TOP_CENTER);
			hbox.setSpacing(PADDING);
			
			FlowPane flow = new FlowPane();
			RenderMatrix render = new RenderMatrix(model, training, false, invertDir);
			render.draw();
			flow.getChildren().addAll(render.getCanvas(), createTable(model, render));
			hbox.getChildren().add(flow);
			
			if (testing.size() > 0)
			{
				flow = new FlowPane();
        		render = new RenderMatrix(model, testing, true, invertDir);
				render.draw();
				flow.getChildren().addAll(render.getCanvas(), createTable(model, render));
        		hbox.getChildren().add(flow);
			}
			
			content.getChildren().add(hbox);
		}
	}
	
	private void updateContent()
	{
		synchronized (mutex)
		{
			if (busy)
			{
				for (FileGroup fg : fileGroups)
				{
					fg.btnChoose.setDisable(true);
					fg.btnDelete.setDisable(true);
					fg.segType.setDisable(true);
					fg.comboField.setDisable(true);
				}
				textFraction.setDisable(true);
    			btnAddFile.setDisable(true);
    			btnLoad.setDisable(true);
    			btnBuild.setDisable(true);
    			btnPredict.setDisable(true);
    			btnSave.setDisable(true);
			}
			else
			{
				for (FileGroup fg : fileGroups)
				{
					fg.btnChoose.setDisable(false);
					fg.btnDelete.setDisable(false);
					fg.segType.setDisable(false);
					fg.comboField.setDisable(fg.radioPred.isSelected());
				}
				textFraction.setDisable(false);
    			btnAddFile.setDisable(false);
    			btnLoad.setDisable(false);
    			btnBuild.setDisable(exec == null || exec.getTraining().size() < 5);
    			btnPredict.setDisable(exec == null || exec.getModel() == null || exec.getPrediction().size() == 0);
    			btnSave.setDisable(exec == null || exec.getModel() == null || exec.getPrediction().size() == 0);
			}
		}
	}
	
	private Node createTable(CompositeModel model, RenderMatrix render)
	{
		int nbins = model.numBins();
		int[] offCounts = render.getOffCounts();
		float[] offPortion = render.getOffPortion();
		float[] offRandom = render.getOffRandom();
		float[] enrichment = render.getEnrichment();
		
		GridPane grid = new GridPane();
		grid.setHgap(PADDING);
		grid.setVgap(PADDING);
		//grid.setGridLinesVisible(true);
		grid.setSnapToPixel(false);
		grid.setPadding(new Insets(2));
		grid.setBorder(new Border(new BorderStroke(Color.BLACK, BorderStrokeStyle.SOLID, CornerRadii.EMPTY, BorderWidths.DEFAULT)));
		//grid.setStyle("-fx-border: 1px solid; -fx-border-color: black;");
		
		for (int n = 1; n <= 4; n++)
		{
			Label label = new Label(n == 1 ? "Count" : n == 2 ? "Portion" : n == 3 ? "Random" : "Enrichment");
			label.setStyle("-fx-font-weight: bold;");
			grid.add(label, n, 0);
		}
		
		for (int n = 0; n < nbins; n++)
		{
			String txt = n == 0 ? "Correct bin" : " + off by " + n;
			grid.add(new Label(txt), 0, n + 1);
			
			grid.add(new Label(String.valueOf(offCounts[n])), 1, n + 1);
			grid.add(new Label(String.format("%.1f%%", 100 * offPortion[n])), 2, n + 1);
			grid.add(new Label(String.format("%.1f%%", 100 * offRandom[n])), 3, n + 1);
			grid.add(new Label(String.format("%.2f", enrichment[n])), 4, n + 1);
		}
		
		for (Node node : grid.getChildren())
		{
			GridPane.setHalignment(node, HPos.CENTER);
			GridPane.setValignment(node, VPos.CENTER);
		}
		
		return grid;
	}
	
	private void changeType(int idx, int type)
	{
		synchronized (mutex)
		{
			if (busy) return;
			session.getFile(idx).type = type;
		}
		updateContent();
	}
	private void changeField(int idx, String value)
	{
		synchronized (mutex)
		{
			if (busy) return;
			session.getFile(idx).field = value;
		}
	}
	private void changeFraction(String value)
	{
		synchronized (mutex)
		{
			if (busy) return;
			try {session.setFraction(Float.parseFloat(value));}
			catch (NumberFormatException ex) {}
		}
	}

	private void actionChooseFile(int idx)
	{
		synchronized (mutex)
		{
			if (busy) return;
		}	
		
		FileGroup fg = fileGroups.get(idx);
	
		boolean isOutput = fg.radioOutput.isSelected();
		String fn = fg.labelFilename.getText();

        FileChooser chooser = new FileChooser();
    	chooser.setTitle(isOutput ? "Save Predictions" : "Open Datafile");
    	if (fn.length() > 0) chooser.setInitialDirectory(new File(fn).getParentFile());
    	
    	File file = isOutput ? chooser.showSaveDialog(stage) : chooser.showOpenDialog(stage);
		if (file == null) return;
		
		synchronized (mutex)
		{
			if (busy) return; // probably not possible
			Session.DataFile df = session.getFile(idx);
			df.filename = file.getPath();
			df.molecules.clear();
			exec = null;
			recreateContent();
		}
	}
	private void actionDeleteFile(int idx)
	{
		synchronized (mutex)
		{
			if (busy) return;
			session.deleteFile(idx);
			recreateContent();
		}
	}
	private void actionAddFile()
	{
		synchronized (mutex)
		{
			if (busy) return;
			Session.DataFile df = new Session.DataFile("", Session.FILE_TRAINING, "");
			session.addFile(df);
			recreateContent();
		}
	}
	private void actionLoad()
	{
		synchronized (mutex)
		{
			if (busy) return;
			busy = true;
		}
		updateContent();
		
		new Thread(() ->
		{
			if (exec == null) exec = new ExecuteSession(session);
		
			for (int n = 0; n < session.numFiles(); n++)
			{
				try {exec.loadFile(n);}
				catch (IOException ex)
				{
					final String fn = session.getFile(n).filename;
			        Platform.runLater(() -> Util.informMessage("Load Failed", "For file [" + fn + "].\nReason: " + ex.getMessage()));
					ex.printStackTrace();
					break;
				}
			}
			exec.partitionMolecules();
		
			synchronized (mutex)
			{
				busy = false;
			}
	        Platform.runLater(() -> recreateContent());
	        
		}).start();
	}
	private void actionBuild()
	{
		synchronized (mutex)
		{
			if (busy) return;
			busy = true;
		}
		updateContent();
		
		new Thread(() ->
		{
			if (exec == null) exec = new ExecuteSession(session);
		
			try {exec.buildModel();}
			catch (Exception ex)
			{
		        Platform.runLater(() -> Util.informMessage("Model Build Failed", "Reason: " + ex.getMessage()));
				ex.printStackTrace();
			}
		
			synchronized (mutex)
			{
				busy = false;
			}
	        Platform.runLater(() -> recreateContent());
	        
		}).start();
	}
	private void actionPredict()
	{
		synchronized (mutex)
		{
			if (busy) return;
			if (exec == null) return;
			List<CompositeModel.Entry> prediction = exec.getPrediction();
			if (prediction.size() == 0)
			{
				Util.informWarning("Prediction", "There are no molecules to predict. Add a new file, and change the type to Prediction.");
				return;
			}
			CompositeModel model = exec.getModel();
			if (model == null) return;
			
    		Stage stage = new Stage();
    		PredictionWindow wnd = new PredictionWindow(stage, model, new ArrayList<CompositeModel.Entry>(prediction));
    		stage.show();
		}
	}
	private void actionSave()
	{
		synchronized (mutex)
		{
			if (busy) return;
		}
	
		if (exec == null) return;
		if (exec.getPrediction().size() == 0)
		{
			Util.informWarning("Save", "There are no molecules to predict. Add a new file, and change the type to Prediction.");
			return;
		}
		CompositeModel model = exec.getModel();
		if (model == null)
		{
			Util.informWarning("Save", "Build a model first.");
			return;
		}
	
		String fn = null, field = "Prediction";
		for (int n = 0; n < session.numFiles(); n++)
		{
			Session.DataFile df = session.getFile(n);
			if (df.type == Session.FILE_OUTPUT)
			{
				fn = df.filename;
				if (df.field.length() > 0) field = df.field;
				break;
			}
		}
		
		if (fn == null)
		{
            FileChooser chooser = new FileChooser();
        	chooser.setTitle("Save Predictions");
        	
        	File file = chooser.showSaveDialog(stage);
    		if (file == null) return;
    		fn = file.getPath();
    		if (!fn.endsWith(".sdf")) fn += ".sdf";
		}
		
		synchronized (mutex)
		{
			busy = true;
		}
		updateContent();

		final String ufn = fn, ufield = field;
		new Thread(() ->
		{
			Exception fail = null;
			try {exec.saveOutput(ufn, ufield);}
			catch (Exception ex) {fail = ex;}
		
			synchronized (mutex)
			{
				busy = false;
			}
			
			Exception ufail = fail;
	        Platform.runLater(() -> 
	        {
	        	recreateContent();

				if (ufail == null)
				{
					Util.informMessage("Saved", "Written predictions to " + ufn);
				}
				else
    			{
    				Util.informMessage("Output Save Failed", "Reason: " + ufail.getMessage());
    				ufail.printStackTrace();
    			}
	        });
	        
		}).start();
	}
}
