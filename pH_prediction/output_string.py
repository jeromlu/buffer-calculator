# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:18:57 2019

@author: jeromlu2
"""


def create_string(buff_name, chemicals, amounts, pH, cond):
    pH = float(pH)

    out_txt = "<h2>Buffer recepie of " + buff_name + "</h2>"
    out_txt = (
        out_txt + '<p><span style="font-size: large;">Mix together:</span></p><ul>'
    )
    for component, txt_amount in zip(chemicals, amounts):
        amount = float(txt_amount.split(" ")[0])
        if amount > 0:
            out_txt = out_txt + "<li>{0} of {1}</li>".format(txt_amount, component)
    out_txt = out_txt + "</ul><p>Predicted pH is&nbsp; {0} ({1} - {2})</p>".format(
        pH, pH - 0.2, pH + 0.2
    )
    out_txt = (
        out_txt + "<p>Measured pH is __________</p>"
        "<p>Predicted conductivity " + cond + "</p><p>Measured"
        "conductivity __________ mS/cm</p><p>&nbsp;</p>"
        "<p>Buffer was filtered through 0.22 &mu;m filter (LOT:____________, KAT:_____________)</p>"
        "<p>Buffer is stored at temperature from 2 to 28 °C</p>"
    )
    return out_txt
