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
    public partial class RegionInformation : Form
    {
        public RegionInformation(string path)
        {
            InitializeComponent();

            pictureBox1.ImageLocation = path;
            pictureBox1.SizeMode = PictureBoxSizeMode.Zoom;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.Close();
        }
    }
}
