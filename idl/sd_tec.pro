pro sd_tec

common rad_data_blk
common radarinfo
common tec_data_blk
common omn_data_blk
common aur_data_blk
common kpi_data_blk
common amp_data_blk



dateSel = [20150409, 20150409];
timeRange = [0700, 0900];


dt_skip_time=5.d ;;; we search data the grd file every 2 min
del_jul=dt_skip_time/1440.d ;;; This is the time step used to read the data --> Selected to be 60 min







;; we need to adjust the xticks with the number of hours plotted...
;; we need to adjust the xticks with the number of hours plotted...
;; we need to adjust the xticks with the number of hours plotted...
num_hours_for_plot = round( ( timeRange[1]/100 ) - ( timeRange[0]/100 ) )


xticks_major_num = 8
xminor_ticks_num = 6
if ( (num_hours_for_plot eq 24.) or (num_hours_for_plot eq 22.)) then begin
        xticks_major_num = num_hours_for_plot/2
        xminor_ticks_num = 6
endif

if ( (num_hours_for_plot eq 23.) ) then begin
        xticks_major_num = 9
        xminor_ticks_num = 6
endif

if ( (num_hours_for_plot ge 13.) and (num_hours_for_plot le 21.) ) then begin
        xticks_major_num = 9
        xminor_ticks_num = 4
endif

if ( (num_hours_for_plot ge 6.) and (num_hours_for_plot le 12.) ) then begin
        xticks_major_num = num_hours_for_plot
        xminor_ticks_num = 6
endif


if ( (num_hours_for_plot ge 3.) and (num_hours_for_plot lt 6.) ) then begin
        xticks_major_num = 2*num_hours_for_plot
        xminor_ticks_num = 6
endif

if ( (num_hours_for_plot gt 1.) and (num_hours_for_plot le 2.) ) then begin
        xticks_major_num = 12
        xminor_ticks_num = 5
endif


if ( (num_hours_for_plot le 1.) ) then begin
        xticks_major_num = 12
        xminor_ticks_num = 5
endif
;; we need to adjust the xticks with the number of hours plotted...
;; we need to adjust the xticks with the number of hours plotted...
;; we need to adjust the xticks with the number of hours plotted...


;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.
;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.
;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.
year_plot=fix(dateSel[0]/1e4)
mndy=double(dateSel[0])-double(year_plot*1e4)
month_plot=fix(mndy/1e2)
day_plot=fix(mndy-month_plot*1e2)
month_list_plot=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

date_in_plot_format=month_list_plot[month_plot-1]+'/'+strtrim(string(day_plot),2)+'/'+strtrim(string(year_plot),2)
;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.
;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.
;;;; Print the date in a proper format on the plot, so get year, month and day from date variable.















rad_fan_ids = [ 208, 209, 206, 207, 204, 205, 33, 32, 40 ]



;; set plot/map parameters
xrangePlot = [-45,45];[-30, 30]
yrangePlot = [-45,45];[-40,-10]
velScale = [-500,500]
tecScale = [0.,20.]
ampScale = [ -0.8, 0.8 ]
cntrMinVal = 0.1
n_levels = 5
coords = "mlt"
omnCharsize = 0.5



tec_read, dateSel
rad_map_read, dateSel
amp_read, dateSel

sfjul,dateSel,timeRange,sjul_search,fjul_search




nele_search=((fjul_search-sjul_search)/del_jul)+1 ;; Num of 2-min times to be searched..
npanels = round((fjul_search-sjul_search)*1440.d/dt_skip_time) + 1


ps_open, '/home/bharatr/Docs/plots/saps-plot-sd-tec-' + strtrim( string(dateSel[0]), 2) + '.ps'

for srch=0,nele_search-1 do begin
        clear_page
        set_format, /sardi
        ;;;Calculate the current jul
        juls_curr=sjul_search+srch*del_jul
        juls_curr_tec = juls_curr + 5.d/1440.d

        sfjul,dateCurrPlot,timeCurrPlot,juls_curr,/jul_to_date
        sfjul,dateCurrTEC,timeCurrTEC,juls_curr_tec,/jul_to_date

        if (timeCurrPlot ne 0) then begin
                timeCurrTEC = timeCurrPlot
        endif

        ;;_position = define_panel(1, 1, 0, 0, aspect=aspect, /bar) 


        rad_load_colortable,/leicester

        mapPos1 = define_panel(1,1.,0,0., aspect=aspect,/bar)

        titleStr = strtrim( string(dateCurrPlot[0]), 2) + "-" + strtrim( string(timeCurrPlot[0]), 2)

        map_plot_panel,date=dateCurrPlot,time=timeCurrPlot,coords=coords,/no_fill,xrange=xrangePlot, $
                yrange=yrangePlot,pos=mapPos1,/isotropic,grid_charsize='0.5',/no_coast,/north, charsize = 1., title = titleStr


        rad_load_colortable, /bw
        ;;plot tec vectors
        tec_median_filter,date=dateCurrTEC,time=timeCurrTEC, threshold=0.10
        overlay_tec_median, date=dateCurrTEC, time=timeCurrTEC, scale=tecScale, coords=coords
        
        rad_load_colortable,/leicester
        ;rad_map_overlay_vectors, date = dateCurrPlot, time=timeCurrPlot, coords = coords, radar_ids=rad_fan_ids, $
        ;                 /no_fov_names, /no_show_Nvc,/no_vector_scale, scale=velScale, symsize=0.25,fixed_color = get_white();215
        
        ;rad_load_colortable, /bw
        ;rad_map_overlay_poes_bnd, dateCurrPlot, timeCurrPlot, coords = coords, $
        ;        fitline_color = get_red(), fitline_style = 3, $
        ;        fitline_thick = 5

        ;rad_map_overlay_poes, dateCurrPlot, timeCurrPlot, coords = coords

        rad_map_overlay_scan, rad_fan_ids, juls_curr, /AJ_filter, coords=coords, scale=velScale, rad_sct_flg_val=2

        ;overlay_coast, coords=coords, jul=juls_curr, /no_fill
        map_overlay_grid, grid_linestyle=0, grid_linethick=2, grid_linecolor=get_gray()

        rad_load_colortable,/bw
        plot_colorbar, 1., 1., 0.25, 0.2, /square,scale=tecScale,legend='Total Electron Content [TECU]', level_format='(f6.2)',param='power',/keep_first_last_label
        rad_load_colortable, /leicester
        plot_colorbar, 1., 1., 0.22, 0.2, /square, scale=velScale, parameter='velocity',/keep_first_last_label, legend='Velocity [m/s]',/left

        

endfor

ps_close, /no_filename


end
