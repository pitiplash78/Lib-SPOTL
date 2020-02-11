using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace SPOTL
{
    public partial class UserControl_OceanLoadingProperties : UserControl
    {
        public delegate void StatusUpdateHandler(object sender, StatusEventArg e);
        public event StatusUpdateHandler OnUpdateStatus;

        public string pathProperties { get; set; }

        public OceanLoadingProperties.GreenFunctionEntry GreensFunction { get; set; }

        public OceanLoadingProperties.OceanModelEntry[] OceanModel
        {
            get
            {
                if (_OceanModel != null)
                    return _OceanModel.ToArray();
                else
                    return null;
            }
            set 
            {
                OceanLoadingProperties.OceanModelEntry[] x = value;
                if (x != null)
                    _OceanModel = new List<OceanLoadingProperties.OceanModelEntry>(x);
            }
        }
        private List<OceanLoadingProperties.OceanModelEntry> _OceanModel = null;

        public OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel
        {
            get
            {
                if (_LocalOceanModel != null)
                    return _LocalOceanModel.ToArray();
                else
                    return null;
            }
            set
            {
                OceanLoadingProperties.LocalOceanModelEntry[] x = value;
                if( x != null)
                    _LocalOceanModel = new List<OceanLoadingProperties.LocalOceanModelEntry>(x);
            }
        }
        private List<OceanLoadingProperties.LocalOceanModelEntry> _LocalOceanModel = null; 

        public bool IsValid { get; private set; }

        #region private variables
        OceanLoadingProperties oclp = null;
        private bool suppress = false;
        #endregion

        public UserControl_OceanLoadingProperties()
        {
            InitializeComponent();
        }

        public UserControl_OceanLoadingProperties(OceanLoadingProperties.GreenFunctionEntry GreensFunction, 
                                               OceanLoadingProperties.OceanModelEntry[] OceanModel,
                                               OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel,
                                               string pathProperties)
        {
            this.pathProperties = pathProperties;
            InitializeComponent();

            setSensor(GreensFunction, OceanModel, LocalOceanModel, pathProperties);
        }

        public void setSensor(OceanLoadingProperties.GreenFunctionEntry GreensFunction,
                              OceanLoadingProperties.OceanModelEntry[] OceanModel,
                              OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel, string pathProperties)
        {
            this.GreensFunction = GreensFunction;
            if (OceanModel != null)
                this._OceanModel = new List<OceanLoadingProperties.OceanModelEntry>(OceanModel);
            else
                this._OceanModel = null;
            if (LocalOceanModel != null)
                this._LocalOceanModel = new List<OceanLoadingProperties.LocalOceanModelEntry>(LocalOceanModel);
            else
                this._LocalOceanModel = null;

            oclp = OceanLoadingProperties.deserialisieren(pathProperties);

            setFormsEarthModel();

            setFormsOceanModel();

            setFormsLocalOceanModel();

            readyForContinue();
        }

        #region Earth model
        private void setFormsEarthModel()
        {
            comboBoxEarthModel.Items.Clear();

            comboBoxEarthModel.Items.Add("<Select Earth model>");
            int selection = 0;

            List<string> gfstrings = new List<string>();
            for (int i = 0; i < oclp.GreenFunctionItems.Length; i++)
            {
                gfstrings.Add(System.IO.Path.GetFileNameWithoutExtension(oclp.GreenFunctionItems[i].FileName));
            }
            gfstrings.Sort();

            for (int i = 0; i < oclp.GreenFunctionItems.Length; i++)
            {
                string tmp = gfstrings[i];
                comboBoxEarthModel.Items.Add(tmp);
                //comboBoxEarthModel.Items.Add(System.IO.Path.GetFileNameWithoutExtension(oclp.GreenFunctionItems[i].FileName));
                if (GreensFunction != null &&
                   System.IO.Path.GetFileNameWithoutExtension(GreensFunction.FileName) == tmp)
                  // System.IO.Path.GetFileNameWithoutExtension(oclp.GreenFunctionItems[i].FileName))
                    selection = i + 1;
            }

            // set Earth Model from SensorEntry
            suppress = true;
            comboBoxEarthModel.SelectedIndex = selection;
            comboBoxEarthModel.Refresh();
            suppress = false;
        }

        private void comboBoxEarthModel_SelectedIndexChanged(object sender, EventArgs e)
        {
            updateSensorEntry();
        }
        #endregion

        #region Ocean model
        private Button[] models = null;
        private CheckBox[] cbWaves = null;

        private void setFormsOceanModel()
        {
            panelOceanModel.Controls.Clear();

            ToolTip toolTip1 = new ToolTip();
            toolTip1.AutoPopDelay = 5000;
            toolTip1.InitialDelay = 1000;
            toolTip1.ReshowDelay = 500;
            toolTip1.ShowAlways = false;

            Button[] waves = new Button[oclp.WaveSequenceOceanModel.Length];

            for (int i = 0; i < oclp.WaveSequenceOceanModel.Length; i++)
            {
                waves[i] = new Button()
                    {
                        Location = new Point(panelOceanModel.Location.X + i * 30 + 100, panelOceanModel.Location.Y - 24),
                        Size = new Size(30, 23),
                        Enabled = false,
                        Text = oclp.WaveSequenceOceanModel[i],
                    };
            }
            models = new Button[oclp.OceanModelItems.Length];
            int numberOfWaves = 0;

            for (int i = 0; i < oclp.OceanModelItems.Length; i++)
            {
                numberOfWaves += oclp.OceanModelItems[i].DarwinSymbol.Length;
                models[i] = new Button()
                  {
                      Location = new Point(0, i * 23),
                      Size = new Size(100, 23),
                      Enabled = true,
                      TextAlign = ContentAlignment.MiddleLeft,
                      Text = oclp.OceanModelItems[i].ModelName,
                      Name = oclp.OceanModelItems[i].ModelName,
                  };
                models[i].Click += new EventHandler(models_Click);
                toolTip1.SetToolTip(models[i], "Description: " + oclp.OceanModelItems[i].Description + Environment.NewLine +
                                               "Year: " + oclp.OceanModelItems[i].Year);
            }

            cbWaves = new CheckBox[numberOfWaves];
            int cbWaveCount = 0;
            for (int i = 0; i < oclp.OceanModelItems.Length; i++)
            {
                for (int j = 0; j < oclp.OceanModelItems[i].DarwinSymbol.Length; j++)
                {
                    for (int k = 0; k < oclp.WaveSequenceOceanModel.Length; k++)
                    {
                        if (oclp.WaveSequenceOceanModel[k].ToLower() == oclp.OceanModelItems[i].DarwinSymbol[j].ToLower()) 
                        {
                            cbWaves[cbWaveCount] = new CheckBox()
                                {
                                    Location = new Point(k * 30 + 105, i * 23 + 2),
                                    Size = new Size(20, 20),
                                    Name = oclp.OceanModelItems[i].ModelName + "_" + oclp.OceanModelItems[i].DarwinSymbol[j],
                                };

                            // set parameter according to the sensor properties
                            if (_OceanModel != null)
                                for (int l = 0; l < _OceanModel.Count; l++)
                                    for (int m = 0; m < _OceanModel[l].DarwinSymbol.Length; m++)
                                        if (_OceanModel[l].ModelName == oclp.OceanModelItems[i].ModelName &&
                                           _OceanModel[l].DarwinSymbol[m] == oclp.OceanModelItems[i].DarwinSymbol[j])
                                            cbWaves[cbWaveCount].Checked = true;

                            cbWaves[cbWaveCount].CheckedChanged += new EventHandler(cbWaves_CheckedChanged);
                            cbWaveCount++;
                        }
                    }
                }
            }

            this.Controls.AddRange(waves);
            this.panelOceanModel.Controls.AddRange(models);
            this.panelOceanModel.Controls.AddRange(cbWaves);
        }

        private void models_Click(object sender, EventArgs e)
        {
            string model = ((Button)sender).Name;

            bool ischecked = false;
            for (int i = 0; i < cbWaves.Length; i++)
            {
                if (cbWaves[i].Name.Split(new char[] { '_' })[0] == model && cbWaves[i].Checked)
                    ischecked = true;

                cbWaves[i].Checked = false;
            }
            bool status = !ischecked;

            suppress = true;
            for (int i = 0; i < cbWaves.Length; i++)
            {
                if (cbWaves[i].Name.Split(new char[] { '_' })[0] == model)
                    cbWaves[i].Checked = status;
                
            }
            suppress = false;

            updateSensorEntry();
        }

        void cbWaves_CheckedChanged(object sender, EventArgs e)
        {
            if (suppress)
                return;

            CheckBox cb = (CheckBox)sender;
            if (cb.Checked)
                for (int i = 0; i < cbWaves.Length; i++)
                    if (cbWaves[i].Name.Split(new char[] { '_' })[1] ==  cb.Name.Split(new char[] { '_' })[1] && 
                        cb.Name != cbWaves[i].Name &&
                        cbWaves[i].Checked)
                        cbWaves[i].Checked = false;

            updateSensorEntry();
        }
        #endregion

        #region Local ocean model
        private Button[] localmodels = null;
        private CheckBox[] cbLocalWaves = null;
      
        private void setFormsLocalOceanModel()
        {
            panelLocalOceanModel.Controls.Clear();

            ToolTip toolTip2 = new ToolTip();
            toolTip2.AutoPopDelay = 5000;
            toolTip2.InitialDelay = 1000;
            toolTip2.ReshowDelay = 500;
            toolTip2.ShowAlways = false;

            Button[] waves = new Button[oclp.WaveSequenceLocalModel.Length];

            for (int i = 0; i < oclp.WaveSequenceLocalModel.Length; i++)
                waves[i] = new Button()
                {
                    Location = new Point(panelLocalOceanModel.Location.X + i * 30 + 120, panelLocalOceanModel.Location.Y - 24),
                    Size = new Size(30, 23),
                    Enabled = false,
                    Text = oclp.WaveSequenceLocalModel[i],
                };

            localmodels = new Button[oclp.LocalOceanModelItem.Length];
            int numberOfWaves = 0;

            for (int i = 0; i < oclp.LocalOceanModelItem.Length; i++)
            {
                numberOfWaves += oclp.LocalOceanModelItem[i].DarwinSymbol.Length;
                localmodels[i] = new Button()
                {
                    Location = new Point(0, i * 23),
                    Size = new Size(120, 23),
                    Enabled = true,
                    TextAlign = ContentAlignment.MiddleLeft,
                    Text =  oclp.LocalOceanModelItem[i].Number + " " + oclp.LocalOceanModelItem[i].ModelName,
                    Name = oclp.LocalOceanModelItem[i].ModelName,
                };
                localmodels[i].Click += new EventHandler(localModels_Click);
                toolTip2.SetToolTip(localmodels[i], "Location: " + oclp.LocalOceanModelItem[i].Location + Environment.NewLine +
                                                    "Year: " + oclp.LocalOceanModelItem[i].Year);
            }

            cbLocalWaves = new CheckBox[numberOfWaves];
            int cbWaveCount = 0;
            for (int i = 0; i < oclp.LocalOceanModelItem.Length; i++)
            {
                for (int j = 0; j < oclp.LocalOceanModelItem[i].DarwinSymbol.Length; j++)
                {
                    for (int k = 0; k < oclp.WaveSequenceLocalModel.Length; k++)
                    {
                        if (oclp.WaveSequenceLocalModel[k].ToLower() == oclp.LocalOceanModelItem[i].DarwinSymbol[j].ToLower()) 
                        {
                            cbLocalWaves[cbWaveCount] = new CheckBox()
                            {
                                Location = new Point(k * 30 + 125, i * 23 + 2),
                                Size = new Size(20, 20),
                                Name = oclp.LocalOceanModelItem[i].ModelName + "_" + oclp.LocalOceanModelItem[i].DarwinSymbol[j],
                            };

                            // set parameter according to the sensor properties
                            if (_LocalOceanModel != null)
                                for (int l = 0; l < _LocalOceanModel.Count; l++)
                                    for (int m = 0; m < _LocalOceanModel[l].DarwinSymbol.Length; m++)
                                        if (_LocalOceanModel[l].ModelName == oclp.LocalOceanModelItem[i].ModelName &&
                                            _LocalOceanModel[l].DarwinSymbol[m] == oclp.LocalOceanModelItem[i].DarwinSymbol[j])
                                        {
                                            cbLocalWaves[cbWaveCount].Checked = true;
                                            break;
                                        }
                            cbLocalWaves[cbWaveCount].CheckedChanged += new EventHandler(cbLocalWaves_CheckedChanged);
                            cbWaveCount++;
                        }
                    }
                }
            }

            // set button style according to the sensor parameter
            string modelTest = cbLocalWaves[0].Name.Split(new char[] { '_' }, StringSplitOptions.RemoveEmptyEntries)[0];
            bool allWavesSelected = true;

            for (int i = 0; i < cbLocalWaves.Length; i++)
            {
                string[] parts = cbLocalWaves[i].Name.Split(new char[] { '_' }, StringSplitOptions.RemoveEmptyEntries);
                string model = parts[0];
                string wave = parts[1];

                if (modelTest == model)
                {
                    if (!cbLocalWaves[i].Checked)
                        allWavesSelected = false;
                }
                else
                {
                    if (allWavesSelected)
                    {
                        for (int j = 0; j < localmodels.Length; j++)
                        {
                            if (localmodels[j].Name == modelTest)
                            {
                                localmodels[j].FlatStyle = FlatStyle.Flat;
                                break;
                            }
                        }
                    }
                    modelTest = model;
                    allWavesSelected = true;
                }
            }
            

            this.Controls.AddRange(waves);
            this.panelLocalOceanModel.Controls.AddRange(localmodels);
            this.panelLocalOceanModel.Controls.AddRange(cbLocalWaves);
        }

        void cbLocalWaves_CheckedChanged(object sender, EventArgs e)
        {
            CheckBox cb = (CheckBox)sender;
            updateSensorEntry();
        }

        private void localModels_Click(object sender, EventArgs e)
        {
            Button bt = (Button)sender;
            if (bt.FlatStyle == FlatStyle.Standard)
                bt.FlatStyle = FlatStyle.Flat;
            else if (bt.FlatStyle == FlatStyle.Flat)
                bt.FlatStyle = FlatStyle.Standard;

            bool check = false;
            if (bt.FlatStyle == FlatStyle.Flat)
                check = true;

            string model = bt.Name;
            suppress = true;

            for (int i = 0; i < cbLocalWaves.Length; i++)
            {
                string wave = cbLocalWaves[i].Name.Split(new char[] { '_' })[0];
                if (wave == model)
                    cbLocalWaves[i].Checked = check;
            }
            suppress = false;

            updateSensorEntry();
        }
        #endregion 

        private void updateSensorEntry()
        {
            if (suppress)
                return;

            GreensFunction = null;
            _OceanModel = null;
            _LocalOceanModel = null;

            // Earth Model
            string selectedModel = null;
            if (comboBoxEarthModel.SelectedItem != null)
            {
                selectedModel = comboBoxEarthModel.SelectedItem.ToString();

                for (int i = 0; i < oclp.GreenFunctionItems.Length; i++)
                {
                    if (selectedModel == System.IO.Path.GetFileNameWithoutExtension(oclp.GreenFunctionItems[i].FileName))
                    {
                        GreensFunction = new OceanLoadingProperties.GreenFunctionEntry(
                            oclp.GreenFunctionItems[i].Description,
                            oclp.GreenFunctionItems[i].FileName,
                            oclp.GreenFunctionItems[i].EarthModel,
                            oclp.GreenFunctionItems[i].Source,
                            oclp.GreenFunctionItems[i].Pattern,
                            oclp.GreenFunctionItems[i].Frame);
                        break;
                    }
                }
            }
            // Ocean Model
            _OceanModel = new List<OceanLoadingProperties.OceanModelEntry>();
            for( int i = 0 ; i < cbWaves.Length ; i++)
            {
                if (cbWaves[i].Checked)
                {
                    string[] parts = cbWaves[i].Name.Split(new char[] { '_' }, StringSplitOptions.RemoveEmptyEntries);
                    string model = parts[0];
                    string wave = parts[1];
                    
                    bool existent = false;
                    for (int j = 0; j < _OceanModel.Count; j++)
                    {
                        if (OceanModel[j].ModelName == model)
                            existent = true;
                    }

                    if (!existent)
                    {
                        for (int k = 0; k < oclp.OceanModelItems.Length; k++)
                        {
                            if (oclp.OceanModelItems[k].ModelName == model)
                            {
                                _OceanModel.Add(new OceanLoadingProperties.OceanModelEntry()
                                    {
                                        ModelName = oclp.OceanModelItems[k].ModelName,
                                        Description = oclp.OceanModelItems[k].Description,
                                        FilePath = oclp.OceanModelItems[k].FilePath,
                                        Year = oclp.OceanModelItems[k].Year,
                                        DarwinSymbol = new string[] { wave },
                                    });
                                break;
                            }
                        }
                    }
                    else
                    {
                        for (int k = 0; k < OceanModel.Length; k++)
                        {
                            if (_OceanModel[k].ModelName == model)
                            {
                                List<string> ds = new List<string>();
                                if(_OceanModel[k].DarwinSymbol != null)
                                    ds = new List<string>(_OceanModel[k].DarwinSymbol);
                                ds.Add(wave);
                                _OceanModel[k].DarwinSymbol = ds.ToArray();
                                break;
                            }
                        }
                    }
                }
            }

            // Local Model
            _LocalOceanModel = new List<OceanLoadingProperties.LocalOceanModelEntry>();
            for (int i = 0; i < cbLocalWaves.Length; i++)
            {
                if (cbLocalWaves[i].Checked)
                {
                    string[] parts = cbLocalWaves[i].Name.Split(new char[] { '_' }, StringSplitOptions.RemoveEmptyEntries);
                    string model = parts[0];
                    string wave = parts[1];

                    bool existent = false;
                    for (int j = 0; j < _LocalOceanModel.Count; j++)
                    {
                        if (_LocalOceanModel[j].ModelName == model)
                            existent = true;
                    }

                    if (!existent)
                    {
                        for (int k = 0; k < oclp.LocalOceanModelItem.Length; k++)
                        {
                            if (oclp.LocalOceanModelItem[k].ModelName == model)
                            {
                                _LocalOceanModel.Add(new OceanLoadingProperties.LocalOceanModelEntry()
                                {
                                    ModelName = oclp.LocalOceanModelItem[k].ModelName,
                                    Description = oclp.LocalOceanModelItem[k].Description,
                                    Location = oclp.LocalOceanModelItem[k].Location,
                                    FilePath = oclp.LocalOceanModelItem[k].FilePath,
                                    Year = oclp.LocalOceanModelItem[k].Year,
                                    DarwinSymbol = new string[] { wave },
                                });
                                break;
                            }
                        }
                    }
                    else
                    {
                        for (int k = 0; k < _LocalOceanModel.Count; k++)
                        {
                            if (_LocalOceanModel[k].ModelName == model)
                            {
                                List<string> ds = new List<string>();
                                if (_LocalOceanModel[k].DarwinSymbol != null)
                                    ds = new List<string>(_LocalOceanModel[k].DarwinSymbol);
                                ds.Add(wave);
                                _LocalOceanModel[k].DarwinSymbol = ds.ToArray();
                                break;
                            }
                        }
                    }
                }
            }
            if (_LocalOceanModel.Count == 0)
                _LocalOceanModel = null;

            readyForContinue();
        }

        private void readyForContinue()
        {
            bool retVal = true;
            if (GreensFunction == null)
                retVal = false;
            if (_OceanModel == null || _OceanModel.Count == 0)
                retVal = false;

            IsValid = retVal;

            // Make sure someone is listening to event
            if (OnUpdateStatus == null)
                return;

            StatusEventArg args = new StatusEventArg(retVal);
            OnUpdateStatus(this, args);
        }

        private void buttonRegionInformation_Click(object sender, EventArgs e)
        {
            RegionInformation ri = new RegionInformation(System.IO.Path.GetDirectoryName(pathProperties) + System.IO.Path.DirectorySeparatorChar + "LocalModels.bmp");
            ri.Show();
        }

        private void buttonGreenFunctionDescription_Click(object sender, EventArgs e)
        {
            GreenFunctionDescription gfd = new GreenFunctionDescription(System.IO.Path.GetDirectoryName(pathProperties) + System.IO.Path.DirectorySeparatorChar + "GreenFunctionDescription.txt");
            gfd.Show();
        }
    }

    public class StatusEventArg : EventArgs
    {
        public bool Status { get; private set; }

        public StatusEventArg(bool Status)
        {
            this.Status = Status;
        }
    }
}
