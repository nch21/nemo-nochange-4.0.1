module cubic_interp
    implicit none
    public :: scurve, itau
contains
    function scurve(x,x0,dx) result(value)
        real, intent(in) :: x,x0,dx
        real :: value
        value = min(1.0, max(0.0, (x-x0)/dx))
        value = value*value*(3.0-2.0*value)

    end function scurve

    function itau(ytau,taud,y) result(tau)
        real, intent(in) :: ytau(:),taud(:)
        real, intent(in) :: y
        real :: tau
        integer :: ks


        do ks = 1, size(ytau)
            if (ytau(ks) > y) then
                exit
            end if
        end do
        ks = ks - 1
        if (ks<1) then
            ks = 1
        end if
        tau = taud(ks) + (taud(ks+1)-taud(ks))* scurve(y,ytau(ks),ytau(ks+1)-ytau(ks))
        

    end function itau
end module cubic_interp