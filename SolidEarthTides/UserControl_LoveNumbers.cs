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
    public partial class UserControl_LoveNumbers : UserControl
    {
        public SPOTL.LoveNumbers LoveNumbers
        {
            get
            {
                return _LoveNumbers;
            }
            set 
            {
                _LoveNumbers = value;
                setValues();
            }
        }
        private SPOTL.LoveNumbers _LoveNumbers = new LoveNumbers();

        public UserControl_LoveNumbers()
        {
            InitializeComponent();
        }

        private void setValues()
        {
            nUD_h1.Value = (decimal)_LoveNumbers.h[0];
            nUD_h2.Value = (decimal)_LoveNumbers.h[1];
            nUD_h3.Value = (decimal)_LoveNumbers.h[2];
            nUD_k1.Value = (decimal)_LoveNumbers.k[0];
            nUD_k2.Value = (decimal)_LoveNumbers.k[1];
            nUD_k3.Value = (decimal)_LoveNumbers.k[2];
            nUD_l1.Value = (decimal)_LoveNumbers.l[0];
            nUD_l2.Value = (decimal)_LoveNumbers.l[1];
            nUD_l3.Value = (decimal)_LoveNumbers.l[2];
        }

        public void SetDefaultValues()
        {
            _LoveNumbers = new LoveNumbers();
            setValues();
            valueChanged(ValueChangedEvent.Default);
        }

        private void nUD_h1_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.h[0] = (double)nUD_h1.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_h2_ValueChanged(object sender, EventArgs e)
        {
           _LoveNumbers.h[1] = (double)nUD_h2.Value;
           valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_h3_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.h[2] = (double)nUD_h3.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_k1_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.k[0] = (double)nUD_k1.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_k2_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.k[1] = (double)nUD_k2.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_k3_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.k[2] = (double)nUD_k3.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_l1_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.l[0] = (double)nUD_l1.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_l2_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.l[1] = (double)nUD_l2.Value;
            valueChanged(ValueChangedEvent.Individual);
        }
        private void nUD_l3_ValueChanged(object sender, EventArgs e)
        {
            _LoveNumbers.l[2] = (double)nUD_l3.Value;
            valueChanged(ValueChangedEvent.Individual);
        }

        public delegate void ValueChangedHandler(object sender, ValueChangedEventArgs e);
        /// <summary>
        /// Event will by thrown do tue a change of the 'Station' object
        /// </summary>
        public event ValueChangedHandler ValueChanged;

        public enum ValueChangedEvent
        {
            None,
            Individual,
            Default,
        }

        private void valueChanged(ValueChangedEvent status)
        {
            // Make sure someone is listening to event
            if (ValueChanged == null) return;

            ValueChangedEventArgs args = new ValueChangedEventArgs(status);
            ValueChanged(this, args);
        }

        public class ValueChangedEventArgs : EventArgs
        {
            public ValueChangedEvent Status { get; private set; }

            public ValueChangedEventArgs(ValueChangedEvent status)
            {
                Status = status;
            }
        }
    }
}
