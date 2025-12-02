module topo_builder
    implicit none
    public :: heaviside, box, cone, clipped_cone, coastal_sprofile, dist_from_line
contains
    function heaviside(x, x0)
        ! Returns 0 for x < x0, 1 for x >= x0
        real, intent(in) :: x(:), x0
        real :: heaviside(size(x))
        heaviside = 0.
        where (x >= x0) heaviside = 1.
    end function heaviside

    function box(x, x0, x1)
        ! Returns 0 for x < x0, 1 or x0 <= x <= x1, 0 for x > x1
        real, intent(in) :: x(:), x0, x1
        real :: box(size(x))
        box = 0.
        where (x >= x0 .and. x <= x1) box = 1.
    end function box

    function cone(x, x0, dx)
        ! Returns 0 for |x-x0| > dx, straight lines peaking at x = x0
        real, intent(in) :: x(:,:), x0, dx
        real :: cone(size(x,1),size(x,2))
        cone = max(0., 1. - abs(x-x0)/dx)
    end function cone

    function clipped_cone(x, x0, dx, clip)
        ! Returns a cone clipped at height 'clip'
        real, intent(in) :: x(:,:), x0, dx, clip
        real :: clipped_cone(size(x,1),size(x,2 ))
        clipped_cone = min(clip, cone(x, x0, dx))
    end function clipped_cone

    function scurve(x, x0, dx)
        ! Returns 0 for x<x0 or x>x+dx, and a cubic in between.
        real, intent(in) :: x(:,:), x0, dx
        real :: scurve(size(x,1),size(x,2))
        real :: s(size(x,1),size(x,2))
        s = min(1., max(0., (x-x0)/dx))
        scurve = (3. - 2.*s)*( s*s )
    end function scurve

    function coastal_sprofile(x, x0, dx, shelf, lf, bf, sf)
        ! A 'coastal profile' with coastal shelf and slope.
        ! Of profile width dx:
        !   - lf is the land fraction (value 0)
        !   - bf is the fraction that is the beach slope.
        !   - sf is the fraction that is the shelf slope.
        ! The remaining fraction is the shelf.
        real, intent(in) :: x(:,:), x0, dx, shelf, lf, bf, sf
        real :: coastal_sprofile(size(x,1),size(x,2))
        real :: s(size(x,1),size(x,2)), sbs(size(x,1),size(x,2)), ssd(size(x,1),size(x,2))
        s = ( x - x0 )/dx
        sbs = s - lf
        ssd = s - (1.-sf)
        coastal_sprofile = shelf * scurve(sbs,0.,bf) + ( 1. - shelf ) * scurve(ssd,0.,sf)
    end function coastal_sprofile

    function dist_from_line(X,x0,Y,y0,y1)
        ! Returns distance from line x=x0 between y=y0 and y=y1
        real, intent(in) :: X(:,:), x0, Y(:,:), y0, y1
        real :: dist_from_line(size(X,1),size(X,2))
        real :: dx(size(X,1),size(X,2)), yr(size(X,1),size(X,2)), dy(size(X,1),size(X,2))
        dx = X - x0
        yr = min( max(Y, y0), y1)
        dy = Y - yr
        dist_from_line = sqrt( dx*dx + dy*dy)
    end function dist_from_line

end module topo_builder


