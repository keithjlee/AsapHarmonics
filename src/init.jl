function sphevaluate(l::Int64, m::Int64; resolution = resolution)
    C = zeros(resolution+1, 2resolution+1)

    C[sph_mode(l,m)] = 1

    sph_evaluate(C)
end

Y = [[sphevaluate(l, m; resolution = resolution) for m = -l:l] for l = 0:l_max]

thetarange = range(0, pi, resolution+1)
phirange = range(0, 2pi, 2resolution+1)

make_xsphere(thetarange::AbstractRange, phirange::AbstractRange) = [cos(phi) * sin(theta) for theta in thetarange, phi in phirange]
make_ysphere(thetarange::AbstractRange, phirange::AbstractRange) = [sin(phi) * sin(theta) for theta in thetarange, phi in phirange]
make_zsphere(thetarange::AbstractRange, phirange::AbstractRange) = [cos(theta) for theta in thetarange, phi in phirange]
export make_xsphere, make_ysphere, make_zsphere


xsphere = [cos(phi) * sin(theta) for theta in thetarange, phi in phirange]
ysphere = [sin(phi) * sin(theta) for theta in thetarange, phi in phirange]
zsphere = [cos(theta) for theta in thetarange, phi in phirange]
export xsphere, ysphere, zsphere

