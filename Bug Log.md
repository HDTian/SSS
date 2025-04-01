### Bug Log

``` r
Warning message:
Failed to fit group -1.
Caused by error in `smooth.construct.cr.smooth.spec()`:
! x has insufficient unique values to support 60 knots: reduce k. 
```

In the PIP polt `geom_smooth` used smoothed curve but sometimes it cannot work; so now use `geom_point()+ geom_line()` and no smooth)

``` r
Error in data.frame(positions = rep(tpoints, times = L_used), PIP = as.vector(t(PIPres)),  : 
  arguments imply differing number of rows: 20, 198
In addition: There were 50 or more warnings (use warnings() to see the first 50)
```

`weight_function( tpoints , selected_dat$Z , selected_dat$X )` does not return the same length vector of tpoints; check if `weight_function` is updated or not

``` r
Error in data.frame(positions = rep(tpoints, times = L_used), PIP = as.vector(t(PIPres)),  : 
  arguments imply differing number of rows: 20, 198
In addition: There were 50 or more warnings (use warnings() to see the first 50)
```

`weight_function( tpoints , selected_dat$Z , selected_dat$X )` does not return the same length vector of tpoints; check if `weight_function` is updated or not
