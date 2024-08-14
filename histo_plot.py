
def histo_plot(XYHist, XZHist, title):
    x = np.linspace(0, 20, 100) #x coordinates (Re)
    y = np.linspace(-30, 30, 100) #y coordinates (Re)
    z = 0 #z coordinates in Re

    [Xgipm,Ygipm,Zgipm] = np.meshgrid(x,y,z,indexing="ij")

    fitting_coeffs = merka05_surface_eq_array_GIPM(XMA_all)

    Xn = Xgipm
    Yn = Ygipm
    Zn = Zgipm
    f = fitting_coeffs[0]*Xn**2 + fitting_coeffs[1]*Yn**2 + fitting_coeffs[2]*Zn**2+ 2*fitting_coeffs[3]*Xn*Yn + 2*fitting_coeffs[4]*Yn*Zn + 2*fitting_coeffs[5]*Xn*Zn + 2*fitting_coeffs[6]*Xn+2*fitting_coeffs[7]*Yn + 2*fitting_coeffs[8]*Zn + fitting_coeffs[9]

    x_1 = np.linspace(0, 20, 100) #x coordinates (Re)
    y_1 = 0 #y coordinates (Re)
    z_1 = np.linspace(-30, 30, 100) #z coordinates in Re

    [Xgipm_1,Ygipm_1,Zgipm_1] = np.meshgrid(x_1,y_1,z_1,indexing="ij")
    Xn_1 = Xgipm_1
    Yn_1 = Ygipm_1
    Zn_1 = Zgipm_1
    f_1 = fitting_coeffs[0]*Xn_1**2 + fitting_coeffs[1]*Yn_1**2 + fitting_coeffs[2]*Zn_1**2+ 2*fitting_coeffs[3]*Xn_1*Yn_1 + 2*fitting_coeffs[4]*Yn_1*Zn_1 + 2*fitting_coeffs[5]*Xn_1*Zn_1 + 2*fitting_coeffs[6]*Xn_1+2*fitting_coeffs[7]*Yn_1 + 2*fitting_coeffs[8]*Zn_1 + fitting_coeffs[9]

    #magnetopause model, D = 2 nPa

    m_1 = 10.22
    m_2 = 1.29
    m_3 = 0.184
    m_4 = 8.14
    m_5 = 6.6
    m_6 = 0.58
    m_7 = -0.007
    m_8 = 0.024

    #use B_z =0
    B_z = 0
    D_p = 2

    alpha = (m_6 + m_7*B_z)*(1 +m_8*(np.log(D_p)))

    tanh_angle = m_3*(B_z+m_4)
    r_0 = (m_1 + m_2*np.tanh(tanh_angle))*(D_p**(-1/m_5))

    pi = np.pi

    theta = np.arange(-pi/2, pi/2, 0.01)

    r_mod = (2/(1+np.cos(theta)))**alpha
    r = r_0*r_mod

    X_shue = r*(np.cos(theta))
    R_shue = r*(np.sin(theta))


    ###################
    fig, ax = plt.subplots()
    subfigs = fig.subfigures(1, 1)
    axsLeft = subfigs.subplots(1, 2, sharey=False)
    subfigs.suptitle(title)

    ax0 = axsLeft[0]

    ax0.contour(Xgipm[:,:,0],Ygipm[:,:,0],f[:,:,0],levels = [0],colors="black",linewidths=1)
    ax0.plot(X_shue, R_shue, linewidth=1, color='k')

    ax0.set_aspect('equal')
    ax0.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    ax0.set_ylabel("$Y_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")


    #want to find three points on the bow shock surface, at y=0 and y=±8 and draw lines from there
    #without exceeding current bounds of plot
    inter_med = fitting_coeffs[6]**2 - (fitting_coeffs[0]*fitting_coeffs[9])
    X_BS_nose = (-fitting_coeffs[6] + np.sqrt(inter_med))/fitting_coeffs[0]

    x_s = X_BS_nose
    y_s = 0
    x_e = 30
    y_e = x_e*(-tan_angle_0)

    #want to also have line for just solar wind flow along y=0

    ax0.hlines(y=0, xmin= 0, xmax=30, linewidth=1, color='k')
    #ax.plot([x_s, x_e], [y_s, y_e], color='k',linewidth=1)
    cmap = matplotlib.colormaps.get_cmap('Blues') 
    cmap.set_bad(color='lightgrey')
    im = ax0.imshow(XYHist, interpolation='nearest', origin='lower', extent=[xedg[0], xedg[-1], yedg[0], yedg[-1]], vmax=10_000, cmap = cmap)
    #ax.set_ylim(-30,30)
    #ax.set_xlim(0,30)
    ax0.invert_xaxis()
    ax0.invert_yaxis()
    fig.colorbar(mappable=im,location='bottom',anchor=(0.5, 0), panchor=(0.5, 0.2), pad=0.1, ax=axsLeft, label='No. of Observations')

    ax1 = axsLeft[1]

    ax1.contour(Xgipm_1[:,0,:],Zgipm_1[:,0,:],f_1[:,0,:],levels = [0],colors="black",linewidths=1)
    ax1.plot(X_shue, R_shue, linewidth=1, color='k')

    ax1.set_aspect('equal')
    ax1.set_xlabel("$X_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")
    ax1.set_ylabel("$Z_\\mathrm{GIPM}$ ($R_\\mathrm{E}$)")


    ax1.hlines(y=0, xmin= 0, xmax=30, linewidth=1, color='k')
    #ax.plot([x_s, x_e], [y_s, y_e], color='k',linewidth=1)
    ax1.imshow(XZHist, interpolation='nearest', origin='lower', extent=[xedg[0], xedg[-1], yedg[0], yedg[-1]], vmax=10_000, cmap = cmap)
    #ax.set_ylim(-30,30)
    #ax.set_xlim(0,30)
    ax1.invert_xaxis()
    ax1.invert_yaxis()
    #fig.show()
