namespace SPOTL
{
    partial class SelectOceanLoadingProperties
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.userControl_OceanLoadingProperties = new SPOTL.UserControl_OceanLoadingProperties();
            this.buttonContinue = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // uControl_OceanLoadingProperties1
            // 
            this.userControl_OceanLoadingProperties.GreensFunction = null;
            this.userControl_OceanLoadingProperties.LocalOceanModel = null;
            this.userControl_OceanLoadingProperties.Location = new System.Drawing.Point(12, 12);
            this.userControl_OceanLoadingProperties.Name = "uControl_OceanLoadingProperties1";
            this.userControl_OceanLoadingProperties.OceanModel = null;
            this.userControl_OceanLoadingProperties.pathProperties = null;
            this.userControl_OceanLoadingProperties.Size = new System.Drawing.Size(519, 363);
            this.userControl_OceanLoadingProperties.TabIndex = 0;
            // 
            // buttonContinue
            // 
            this.buttonContinue.Location = new System.Drawing.Point(456, 394);
            this.buttonContinue.Name = "buttonContinue";
            this.buttonContinue.Size = new System.Drawing.Size(75, 23);
            this.buttonContinue.TabIndex = 1;
            this.buttonContinue.Text = "Continue";
            this.buttonContinue.UseVisualStyleBackColor = true;
            this.buttonContinue.Click += new System.EventHandler(this.buttonContinue_Click);
            // 
            // label1
            // 
            this.label1.ForeColor = System.Drawing.Color.Red;
            this.label1.Location = new System.Drawing.Point(12, 378);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(519, 13);
            this.label1.TabIndex = 2;
            // 
            // SelectOceanLoadingProperties
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(548, 428);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.buttonContinue);
            this.Controls.Add(this.userControl_OceanLoadingProperties);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedToolWindow;
            this.Name = "SelectOceanLoadingProperties";
            this.Text = "Select Ocean Loading Properties";
            this.ResumeLayout(false);

        }

        #endregion

        private SPOTL.UserControl_OceanLoadingProperties userControl_OceanLoadingProperties;
        private System.Windows.Forms.Button buttonContinue;
        private System.Windows.Forms.Label label1;
    }
}