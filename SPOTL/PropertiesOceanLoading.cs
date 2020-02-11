using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;

namespace SPOTL
{
    public partial class PropertiesOceanLoading : Form
    {
        public string oclEarthModel = null;
        public string oclTideModel = null;
        private OceanLoadingProperties oclP;

        public PropertiesOceanLoading(OceanLoadingProperties oclP,
                                      string oclEarthModel,
                                      string oclTideModel)
        {
            this.oclP = oclP;
            this.oclEarthModel = oclEarthModel;
            this.oclTideModel = oclTideModel;

            InitializeComponent();

            int tmp = 0;
            int c = 0;
            foreach (OceanLoadingProperties.OceanModelEntry tmname in oclP.OceanModelItems)
            {
                comboBoxOceanLoad.Items.Add(tmname.ModelName);
                if (tmname.FilePath == oclTideModel)
                    tmp = c;
                c++;
            }
            comboBoxOceanLoad.SelectedIndex = tmp;
            tmp = 0;
            c = 0;

            foreach (OceanLoadingProperties.GreenFunctionEntry gfname in oclP.GreenFunctionItems)
            {
                comboBoxGreenFuntion.Items.Add(gfname.Description);
                if (gfname.FileName == oclEarthModel)
                    tmp = c;
                c++;
            }
            comboBoxGreenFuntion.SelectedIndex = tmp;
        }

        private void buttonOK_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.OK;
            this.Close();
        }

        private void buttonCancel_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.OK;
            this.Close();
        }

        private void comboBoxOceanLoad_SelectedIndexChanged(object sender, EventArgs e)
        {
            oclTideModel = oclP.OceanModelItems[comboBoxOceanLoad.SelectedIndex].ModelName;
        }

        private void comboBoxGreenFuntion_SelectedIndexChanged(object sender, EventArgs e)
        {
            oclEarthModel = oclP.GreenFunctionItems[comboBoxGreenFuntion.SelectedIndex].Description;
        }
    }
}
