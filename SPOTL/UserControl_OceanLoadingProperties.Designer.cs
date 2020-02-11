namespace SPOTL
{
    partial class UserControl_OceanLoadingProperties
    {
        /// <summary> 
        /// Erforderliche Designervariable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary> 
        /// Verwendete Ressourcen bereinigen.
        /// </summary>
        /// <param name="disposing">True, wenn verwaltete Ressourcen gelöscht werden sollen; andernfalls False.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Vom Komponenten-Designer generierter Code

        /// <summary> 
        /// Erforderliche Methode für die Designerunterstützung. 
        /// Der Inhalt der Methode darf nicht mit dem Code-Editor geändert werden.
        /// </summary>
        private void InitializeComponent()
        {
            this.comboBoxEarthModel = new System.Windows.Forms.ComboBox();
            this.panelOceanModel = new System.Windows.Forms.Panel();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.panelLocalOceanModel = new System.Windows.Forms.Panel();
            this.buttonRegionInformation = new System.Windows.Forms.Button();
            this.buttonGreenFunctionDescription = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // comboBoxEarthModel
            // 
            this.comboBoxEarthModel.FormattingEnabled = true;
            this.comboBoxEarthModel.Location = new System.Drawing.Point(70, 0);
            this.comboBoxEarthModel.Name = "comboBoxEarthModel";
            this.comboBoxEarthModel.Size = new System.Drawing.Size(377, 21);
            this.comboBoxEarthModel.TabIndex = 0;
            this.comboBoxEarthModel.SelectedIndexChanged += new System.EventHandler(this.comboBoxEarthModel_SelectedIndexChanged);
            // 
            // panelOceanModel
            // 
            this.panelOceanModel.AutoScroll = true;
            this.panelOceanModel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panelOceanModel.Location = new System.Drawing.Point(1, 49);
            this.panelOceanModel.Name = "panelOceanModel";
            this.panelOceanModel.Size = new System.Drawing.Size(516, 141);
            this.panelOceanModel.TabIndex = 1;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(-2, 3);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(66, 13);
            this.label1.TabIndex = 2;
            this.label1.Text = "Earth model:";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(-2, 31);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(73, 13);
            this.label2.TabIndex = 3;
            this.label2.Text = "Ocean model:";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(-2, 202);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(67, 13);
            this.label3.TabIndex = 5;
            this.label3.Text = "Local model:";
            // 
            // panelLocalOceanModel
            // 
            this.panelLocalOceanModel.AutoScroll = true;
            this.panelLocalOceanModel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panelLocalOceanModel.Location = new System.Drawing.Point(1, 220);
            this.panelLocalOceanModel.Name = "panelLocalOceanModel";
            this.panelLocalOceanModel.Size = new System.Drawing.Size(516, 141);
            this.panelLocalOceanModel.TabIndex = 4;
            // 
            // buttonRegionInformation
            // 
            this.buttonRegionInformation.Location = new System.Drawing.Point(450, 196);
            this.buttonRegionInformation.Name = "buttonRegionInformation";
            this.buttonRegionInformation.Size = new System.Drawing.Size(66, 23);
            this.buttonRegionInformation.TabIndex = 6;
            this.buttonRegionInformation.Text = "Region ...";
            this.buttonRegionInformation.UseVisualStyleBackColor = true;
            this.buttonRegionInformation.Click += new System.EventHandler(this.buttonRegionInformation_Click);
            // 
            // buttonGreenFunctionDescription
            // 
            this.buttonGreenFunctionDescription.Location = new System.Drawing.Point(450, 0);
            this.buttonGreenFunctionDescription.Name = "buttonGreenFunctionDescription";
            this.buttonGreenFunctionDescription.Size = new System.Drawing.Size(66, 23);
            this.buttonGreenFunctionDescription.TabIndex = 7;
            this.buttonGreenFunctionDescription.Text = "More ...";
            this.buttonGreenFunctionDescription.UseVisualStyleBackColor = true;
            this.buttonGreenFunctionDescription.Click += new System.EventHandler(this.buttonGreenFunctionDescription_Click);
            // 
            // uControl_OceanLoadingProperties
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.buttonGreenFunctionDescription);
            this.Controls.Add(this.buttonRegionInformation);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.panelLocalOceanModel);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.panelOceanModel);
            this.Controls.Add(this.comboBoxEarthModel);
            this.Name = "uControl_OceanLoadingProperties";
            this.Size = new System.Drawing.Size(519, 363);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.ComboBox comboBoxEarthModel;
        private System.Windows.Forms.Panel panelOceanModel;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Panel panelLocalOceanModel;
        private System.Windows.Forms.Button buttonRegionInformation;
        private System.Windows.Forms.Button buttonGreenFunctionDescription;
    }
}
