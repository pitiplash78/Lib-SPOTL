using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace SPOTL
{
    public partial class SelectOceanLoadingProperties : Form
    {

        public OceanLoadingProperties.GreenFunctionEntry GreensFunction {get; private set;}
        public OceanLoadingProperties.OceanModelEntry[] OceanModel  {get; private set;}
        public OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel { get; private set; }
       // public string pathProperties { get; private set; }

        public SelectOceanLoadingProperties(OceanLoadingProperties.GreenFunctionEntry GreensFunction,
                                               OceanLoadingProperties.OceanModelEntry[] OceanModel,
                                               OceanLoadingProperties.LocalOceanModelEntry[] LocalOceanModel,
                                               string pathProperties)
        {
            this.GreensFunction = GreensFunction;
            this.OceanModel = OceanModel;
            this.LocalOceanModel = LocalOceanModel;

            InitializeComponent();

            userControl_OceanLoadingProperties.setSensor(GreensFunction, OceanModel, LocalOceanModel, pathProperties);

            userControl_OceanLoadingProperties.OnUpdateStatus += new SPOTL.UserControl_OceanLoadingProperties.StatusUpdateHandler(uControl_OceanLoadingProperties1_OnUpdateStatus);
        }

        void uControl_OceanLoadingProperties1_OnUpdateStatus(object sender, SPOTL.StatusEventArg e)
        {
            continueCheck();
        }

        private void buttonContinue_Click(object sender, EventArgs e)
        {
            GreensFunction = userControl_OceanLoadingProperties.GreensFunction;
            OceanModel = userControl_OceanLoadingProperties.OceanModel;
            LocalOceanModel = userControl_OceanLoadingProperties.LocalOceanModel;
            DialogResult = System.Windows.Forms.DialogResult.OK;
            this.Close();
        }

        private void continueCheck()
        {
            if (!userControl_OceanLoadingProperties.IsValid)
            {
                label1.Text = "Finish the properties and continue!";
                buttonContinue.Enabled = false;
            }
            else
            {
                label1.Text = "";
                buttonContinue.Enabled = true;
            }
        }
    }
}
