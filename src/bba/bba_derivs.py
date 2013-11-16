from math import sin, cos

def PhiDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)) + (-zo_c + zo_pt)*(-sin(ka)*sin(om)*sin(ph) - cos(ka)*cos(ph)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*(-(-xo_c + xo_pt)*cos(om)*cos(ph) + (-zo_c + zo_pt)*sin(ph)*cos(om))*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2

def PhiDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)) + (-zo_c + zo_pt)*(sin(ka)*cos(ph) - sin(om)*sin(ph)*cos(ka)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*(-(-xo_c + xo_pt)*cos(om)*cos(ph) + (-zo_c + zo_pt)*sin(ph)*cos(om))*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2

def OmegaDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*((-xo_c + xo_pt)*sin(om)*sin(ph) + (-yo_c + yo_pt)*cos(om) + (-zo_c + zo_pt)*sin(om)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 - f*((-xo_c + xo_pt)*sin(ka)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(ka)*sin(om) + (-zo_c + zo_pt)*sin(ka)*cos(om)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))
def OmegaDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*((-xo_c + xo_pt)*sin(om)*sin(ph) + (-yo_c + yo_pt)*cos(om) + (-zo_c + zo_pt)*sin(om)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 - f*((-xo_c + xo_pt)*sin(ph)*cos(ka)*cos(om) - (-yo_c + yo_pt)*sin(om)*cos(ka) + (-zo_c + zo_pt)*cos(ka)*cos(om)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))

def KappaDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))
def KappaDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(-sin(ka)*sin(om)*sin(ph) - cos(ka)*cos(ph)) - (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(-sin(ka)*sin(om)*cos(ph) + sin(ph)*cos(ka)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))


def XEODerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(-sin(ka)*sin(om)*sin(ph) - cos(ka)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*sin(ph)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2
def XEODerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(sin(ka)*cos(ph) - sin(om)*sin(ph)*cos(ka))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*sin(ph)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2

def YEODerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*sin(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 + f*sin(ka)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))
def YEODerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*sin(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 + f*cos(ka)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))

def ZEODerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(-sin(ka)*sin(om)*cos(ph) + sin(ph)*cos(ka))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*cos(om)*cos(ph)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2
def ZEODerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(-sin(ka)*sin(ph) - sin(om)*cos(ka)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) - f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*cos(om)*cos(ph)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2


def XPtDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) + f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*sin(ph)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2
def XPtDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) + f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*sin(ph)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2

def YPtDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*sin(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 - f*sin(ka)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))
def YPtDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*sin(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2 - f*cos(ka)*cos(om)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))

def ZPtDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) + f*((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))*cos(om)*cos(ph)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2
def ZPtDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -f*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph)) + f*((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))*cos(om)*cos(ph)/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))**2


def XCenterDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return 1
def XCenterDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return 0

def YCenterDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return 0
def YCenterDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return 1


def FocalDerX(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -((-xo_c + xo_pt)*(sin(ka)*sin(om)*sin(ph) + cos(ka)*cos(ph)) + (-yo_c + yo_pt)*sin(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(om)*cos(ph) - sin(ph)*cos(ka)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))
def FocalDerY(xi_c, yi_c, f, xo_c, yo_c, zo_c, ph, om, ka, xo_pt, yo_pt, zo_pt, **kwargs):
    return -((-xo_c + xo_pt)*(-sin(ka)*cos(ph) + sin(om)*sin(ph)*cos(ka)) + (-yo_c + yo_pt)*cos(ka)*cos(om) + (-zo_c + zo_pt)*(sin(ka)*sin(ph) + sin(om)*cos(ka)*cos(ph)))/((-xo_c + xo_pt)*sin(ph)*cos(om) - (-yo_c + yo_pt)*sin(om) + (-zo_c + zo_pt)*cos(om)*cos(ph))

x_derivs = { "ph" : PhiDerX,
             "om" : OmegaDerX,
             "ka" : KappaDerX,
             "xo_c" : XEODerX,
             "yo_c" : YEODerX,
             "zo_c" : ZEODerX,
             "xo_pt" : XPtDerX,
             "yo_pt" : YPtDerX,
             "zo_pt" : ZPtDerX,
             "xi_c" : XCenterDerX,
             "yi_c" : YCenterDerX,
             "f" : FocalDerX,
           }

y_derivs = { "ph" : PhiDerY,
             "om" : OmegaDerY,
             "ka" : KappaDerY,
             "xo_c" : XEODerY,
             "yo_c" : YEODerY,
             "zo_c" : ZEODerY,
             "xo_pt" : XPtDerY,
             "yo_pt" : YPtDerY,
             "zo_pt" : ZPtDerY,
             "xi_c" : XCenterDerY,
             "yi_c" : YCenterDerY,
             "f" : FocalDerY,
           }
