using System;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Windows.Forms;

namespace SPOTL
{
    public partial class GreenFunctionDescription : Form
    {
        public GreenFunctionDescription(string path)
        {
            InitializeComponent();

            richTextBox1.Text = "";

            if (System.IO.File.Exists(path))
            {
                System.IO.StreamReader re = new System.IO.StreamReader(new System.IO.FileStream(path, System.IO.FileMode.Open));
                richTextBox1.Text = re.ReadToEnd();
                re.Close();
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.Close();
        }
    }
}
