# cd("C:\\Users\\pi96doc\\Documents\\Programming\\Julia\\OrthoView.jl\\")
# ]activate .
## interactive example:
module OrthoView
export ortho!, ortho, test_ortho_view

using GLMakie
using GLMakie.FileIO
using Printf
using ColorTypes  # for RGB type
using ColorSchemes
using PerceptualColourMaps

function get_crosshair_xs(px, sx, gap)
    xs = [0f0,px-gap,px+gap,sx,px,px,px,px]
    return xs
end
function get_crosshair_ys(py, sy, gap)
    gap = Float32(to_value(gap))
    ys = [py,py,py,py,0f0,py-gap,py+gap,sy]
    return ys
end

function crosshair(ax, sl_x,sl_y, sz; markersize=0.7, color=:red)
    if isa(markersize, Slider)
        cross_xs = @lift(get_crosshair_xs($(sl_x.value), sz[1], Float32($(markersize.value))))
        cross_ys = @lift(get_crosshair_ys($(sl_y.value), sz[2], Float32($(markersize.value))))
    else
        gap = Float32(to_value(markersize))
        cross_xs = @lift(get_crosshair_xs($(sl_x.value), sz[1], gap))
        cross_ys = @lift(get_crosshair_ys($(sl_y.value), sz[2], gap))
    end
    linesegments!(ax, cross_xs, cross_ys, color=color, linewidth=3)
end

function fit_into(start, stop, sz, do_limit=false)
    if start<0
        stop -= start; start=0;
    end
    if start>sz+1
        stop -= start-(sz+1); start = sz+1;
    end
    if stop<0
        start -= stop; stop=0;
    end
    if stop>sz+1
        start -= stop-(sz+1); stop=sz+1;
    end
    if abs(stop-start)+ 1 >sz
        delta = sign(stop-start)*(abs(stop-start) - (sz))
        start += delta/2 
        stop += delta/2 
    end
    if do_limit
        return clamp(start,-0.5,sz), clamp(stop,-0.5,sz)
    else
        return start, stop
    end
end

function do_zoom(ax_xy, ctr, sz, abs_zoom, aspects=(1,1,1), indices=(1,2))
    if isnothing(ax_xy)
        return
    end
    ix,iy = abs.(indices)
    sz2d = (sz[ix],sz[iy])
    # lx,ly = to_value(ax_xy.limits)
    sx,sy = ax_xy.layoutobservables.computedbbox.val.widths  # The size in pixels that is available
    #if isnothing(lx)
    #    return Consume(true)
    #end
    c = (to_value(ctr) .- 0.5) # .* aspects
    c = (c[ix], c[iy])
    half_width = sign(indices[1])*abs_zoom .* sx / aspects[ix] / 2  
    nx_start, nx_stop = fit_into(c[1] - half_width, c[1] + half_width, sz2d[1], false)
    xlims!(ax_xy,nx_start,nx_stop)
    half_width = sign(indices[2])*abs_zoom .* sy / aspects[iy] / 2  
    ny_start, ny_stop = fit_into(c[2] - half_width, c[2] + half_width, sz2d[2], false)
    ylims!(ax_xy,ny_start,ny_stop)
end

function get_max_zoom(sz, ax, ax_zy, ax_xz)
    abs_zoom_x = sz[1] / ax.layoutobservables.computedbbox.val.widths[1]
    abs_zoom_y = sz[2] / ax.layoutobservables.computedbbox.val.widths[2]
    if length(sz) > 2 && !isnothing(ax_zy) && !isnothing(ax_xz)
        abs_zoom_z = sz[3] / ax_zy.layoutobservables.computedbbox.val.widths[1]
        abs_zoom_z2 = sz[3] / ax_xz.layoutobservables.computedbbox.val.widths[2]
        return max(abs_zoom_x,abs_zoom_y, abs_zoom_z, abs_zoom_z2)
    else
        return max(abs_zoom_x,abs_zoom_y)
    end
end

"""
    zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)

    zooms in or out equal amount with all three connected panels
"""
function zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
    if !isnothing(event)
        # zoom-in or zoom-out
        zoom_dir = to_value(event.y)
        zoom_max = get_max_zoom(sz .* aspects, ax, ax_zy, ax_xz)
        abs_zoom = min(abs_zoom * 1.1 ^ -zoom_dir, zoom_max)
    end
    do_zoom(ax, ctr, sz, abs_zoom, aspects, (1,-2))
    do_zoom(ax_xz, ctr, sz, abs_zoom, aspects, (1,-3))
    do_zoom(ax_zy, ctr, sz, abs_zoom, aspects, (3,-2))
    return abs_zoom
end

function on_color_triangle(cb)
    scene = cb.parent.scene
    mouse = scene.events.mouseposition[]
    c_start = cb.layoutobservables.suggestedbbox.val.origin
    c_stop = cb.layoutobservables.suggestedbbox.val.origin .+ cb.layoutobservables.suggestedbbox.val.widths
    d_top = abs(mouse[2] .- c_start[2])
    d_buttom = abs(mouse[2] .- c_stop[2])
    if !(all(mouse.>c_start) && all(mouse .<c_stop))
        return 0;
    elseif (d_top < 20)
        return -1;
    elseif (d_buttom < 20)
        return 1;
    else
        return 0;
    end
end

function in_colorbar(cb)
    scene = cb.parent.scene
    mouse = scene.events.mouseposition[]
    c_start = cb.layoutobservables.suggestedbbox.val.origin
    c_stop = cb.layoutobservables.suggestedbbox.val.origin .+ cb.layoutobservables.suggestedbbox.val.widths
    return (all(mouse.>c_start) && all(mouse .<c_stop))
end

function register_colorbar_scroll!(cb, colorrange, default_range=(0.0,1.0), high_col=:red, low_col=:blue)
    scene = cb.parent.scene
    on(events(scene).scroll, priority = 2) do event
        if ! (in_colorbar(cb))
            return Consume(false)
        end 
        c_min,c_max = to_value(colorrange)
        zoom_dir = to_value(event[2]) # event.y
        c_rng = (c_max - c_min) * 1.1 ^ -zoom_dir
        if ispressed(scene,Mouse.middle) || ispressed(scene, Keyboard.left_shift) || ispressed(scene, Keyboard.right_shift) || (on_color_triangle(cb) == -1)
            c_min = c_max - c_rng
        else
            c_max = c_min + c_rng
        end
        colorrange[] = (c_min,c_max)
        return Consume(true) # for now prevent the zoom
    end
    on(events(scene).mousebutton, priority = 2) do event
        if ! (in_colorbar(cb))
            return Consume(false)
        end 
        if event.button == Mouse.left && event.action == Mouse.release
            if on_color_triangle(cb) == 1
                if isnothing(to_value(cb.highclip))
                    cb.highclip[] = high_col
                else
                    cb.highclip[] = nothing
                end
                return Consume(true) # for now prevent the zoom
            end
            if on_color_triangle(cb) == -1
                if isnothing(to_value(cb.lowclip))
                    cb.lowclip[] = low_col
                else
                    cb.lowclip[] = nothing
                end
                return Consume(true) # for now prevent the zoom
            end
            if ispressed(scene, Keyboard.left_shift) || ispressed(scene, Keyboard.right_shift)
                max_r = max(default_range...)
                colorrange[] = (-max_r,max_r)
            elseif ispressed(scene, Keyboard.left_control) || ispressed(scene, Keyboard.left_control)
                colorrange[] = (0.0,default_range[2])
            else 
                colorrange[] = default_range
            end
        end
        return Consume(true) # for now prevent the zoom
    end
end

function register_panel_zoom_link!(ax, ctr, sz, ax_xz=nothing, ax_zy=nothing; aspects=ones(length(sz)))
    abs_zoom = get_max_zoom(sz.*aspects, ax, ax_zy, ax_xz)
    zoom_all(abs_zoom, nothing, ax, ax_xz, ax_zy, aspects, ctr, sz)
    register_interaction!(ax, :my_scroll1_interaction) do event::ScrollEvent, axis
        abs_zoom = zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
        return Consume(true) # for now prevent the zoom
    end
    if !isnothing(ax_xz)
        register_interaction!(ax_xz, :my_scroll2_interaction) do event::ScrollEvent, axis
            abs_zoom = zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
            return Consume(true) # for now prevent the zoom
        end
    end
    if !isnothing(ax_zy)
        register_interaction!(ax_zy, :my_scroll2_interaction) do event::ScrollEvent, axis
            abs_zoom = zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
            return Consume(true) # for now prevent the zoom
        end
    end
end

function get_text(name, pos, data, aspects)
    d = data[pos...];
    sz = size(data)
    rel_pos = pos .- (sz .÷2 .+1)
    str = "$(name): "
    if isa(d, Complex)
        str *= "$(@sprintf("%.3f + %.3fi",real(d),imag(d)))\n"
    else
        if isa(d, Real)
            str *= "$(@sprintf("%.3f",d)),\n"
        else
            if isa(d, Tuple)
                d = round.(data[pos...]; digits=3)
                str *= "$(d),\n"
            elseif isa(d,RGB)
                str *= "RGB: $(@sprintf("%.2f",d.r)), $(@sprintf("%.2f",d.g)), $(@sprintf("%.2f",d.b)),\n"
            else
                str *= "$(d),\n"
            end
        end
    end
    str *= "pos: $(pos)\nctr: $(rel_pos)\nsc. ctr: $(rel_pos .* aspects)\nsz:$(size(data))"

end

"""
    register_panel_interactions!(ax, sl_x, sl_y, sl_z, ref_ax; key_buffer="")

    registers the individual interactions which are identical for each panel.
"""
function register_panel_interactions!(ax, sl_x, sl_y, sl_z, ref_ax; key_buffer="")
    sz = (max(to_value(sl_x.range)[1],to_value(sl_x.range)[end]), max(to_value(sl_y.range)[1],to_value(sl_y.range)[end])) # (to_value(ax.limits)[1][2] , to_value(ax.limits)[2][2])
    ctr = sz .÷ 2 .+ 1
    # pos = Observable([Point2f0(0)])
    deregister_interaction!(ax, :rectanglezoom)

    register_interaction!(ax, :my_mouse_interaction) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick || event.type === MouseEventTypes.leftdragstart || 
            event.type === MouseEventTypes.leftdragstop || event.type === MouseEventTypes.leftdrag || 
            event.type === MouseEventTypes.middledragstart || event.type === MouseEventTypes.middledragstop || event.type === MouseEventTypes.middledrag
            set_close_to!(sl_x, event.data[1] )
            set_close_to!(sl_y, event.data[2] )
            ref_ax[1] = ax
        end
    end

    register_interaction!(ax, :my_key_interaction) do event::KeysEvent, axis
        if ref_ax[1] != ax
            return
        end
        if Keyboard.left in event.keys 
            set_close_to!(sl_x, to_value(sl_x.value) - 1f0) 
        end
        if Keyboard.right in event.keys 
            set_close_to!(sl_x, to_value(sl_x.value) + 1f0) 
        end
        if Keyboard.up in event.keys 
            set_close_to!(sl_y, to_value(sl_y.value) - 1f0) 
        end
        if Keyboard.down in event.keys 
            set_close_to!(sl_y, to_value(sl_y.value) + 1f0) 
        end
        if Keyboard.home in event.keys 
            set_close_to!(sl_x, ctr[1]) 
            set_close_to!(sl_y, ctr[2]) 
        end
        if !isnothing(sl_z)
            if Keyboard.page_up in event.keys 
                set_close_to!(sl_z, to_value(sl_z.value) - 1f0) 
            end
            if Keyboard.page_down in event.keys
                set_close_to!(sl_z, to_value(sl_z.value) + 1f0) 
            end
            if Keyboard.enter in event.keys
                try
                    set_close_to!(sl_z, parse(Float32,key_buffer)) 
                    println(key_buffer)
                catch ev
                end
                key_buffer = ""                
            end
            Keyboard._0 in event.keys && (key_buffer *= '0')
            Keyboard._1 in event.keys && (key_buffer *= '1')
            Keyboard._2 in event.keys && (key_buffer *= '2')
            Keyboard._3 in event.keys && (key_buffer *= '3')
            Keyboard._4 in event.keys && (key_buffer *= '4')
            Keyboard._5 in event.keys && (key_buffer *= '5')
            Keyboard._6 in event.keys && (key_buffer *= '6')
            Keyboard._7 in event.keys && (key_buffer *= '7')
            Keyboard._8 in event.keys && (key_buffer *= '8')
            Keyboard._9 in event.keys && (key_buffer *= '9')
        end
    end
end

function get_pos(pos)
        Tuple(round(Int,p) for p in pos)
end

function apply_min_max_gamma(min, max, gamma)
    (dat) -> ((dat - min) / (max-min))^gamma
end

function get_slice(data, pos, dims=(1,2)) # , min_max_val=extrema(abs2.(data))
    idx = Tuple((d in dims) ? Colon() : pos[d] for d = 1:length(pos))

    if dims[1] > dims[2]
        res = transpose(@view data[idx...])
    else
        res = @view data[idx...]
    end
    if eltype(res)<: Complex
        return abs.(res)
    else
        return res
    end
    # res = map(apply_min_max_gamma(min_max_val[1], min_max_val[2], 1.0f0), res)
end

function obs_from_sliders(sls)
    sz = length(sls)
    if sz == 0
        return nothing
    elseif sz == 1
        position = @lift(get_pos(($(sls[1].value))))
    elseif sz == 2
        position = @lift(get_pos(($(sls[1].value), $(sls[2].value))))
    elseif sz == 3
        position = @lift(get_pos(($(sls[1].value), $(sls[2].value), $(sls[3].value))))
    elseif sz == 4
        position = @lift(get_pos(($(sls[1].value), $(sls[2].value), $(sls[3].value), $(sls[4].value))))
    elseif sz == 5
        position = @lift(get_pos(($(sls[1].value), $(sls[2].value), $(sls[3].value), $(sls[4].value), $(sls[5].value))))
    else
        error("a maximum of 5 extra dimensions is allowed")
    end
    position
end

function apply_gamma(rgba, gamma) 
    ColorTypes.RGBA(rgba.r^gamma, rgba.g^gamma, rgba.b^gamma, rgba.alpha)
end

function colormap_line(subfig; textsize=textsize)
    col_options = ["L1","L2","L3","L4","C1","C2","C3","R1","R2","R3"]
    N_cmap = 2048
    my_rawmap = Observable(RGBA{Float32}.(cmap("L1", N=N_cmap)))

    lp_menu = subfig[] = GridLayout()
    ## fill the lower menu bar with content:
    # label=["Color"], 
    menu = Menu(lp_menu[1,1], options = col_options, prompt="Colormap", textsize=textsize) 
    on(menu.selection) do s
        my_rawmap[] = RGBA{Float32}.(cmap(s, N=N_cmap))
    end
    # Gamma slider
    Label(lp_menu[1,2],"Gamma:", halign = :left, valign = :center, textsize=textsize) # 5 Mb
    # , label="Gamma"
    # , align=Inside()
    sg = Slider(lp_menu[1,3], range = -2:0.1:2, horizontal = true, startvalue = 0.0)
    gamma_txt = @lift( "$(@sprintf("%.2f",10.0^$(sg.value)))")
    Label(lp_menu[1,4],gamma_txt, halign = :left, valign = :center, textsize=textsize)
    gamma = @lift(Float32(10f0 ^ $(sg.value)))
    return @lift(apply_gamma.(to_value($(my_rawmap)), to_value($(gamma))))
end

# Memory usage:  Menu: 380 Mb
# Slider, laber, lift: each 5Mb
# hidedecorations: 132 Mb
# axis: 27 Mb
# heatmap: 2 Mb
# crosshair: 0.26 Mb
# Colorbar: 48 Mb
# colsize!: 8.6 Mb
# rowsize!: 8.6 Mb
# xlims!: 1 Mb
# ylims!: 2 Mb


function ortho!(fig, myim; title = "Image", color=:red, markersize = 40.0, aspects=ones(ndims(myim)), colorbar=true, textsize=24)
    sz = size(myim)

    # arguments to hide axes and disable zooming
    hideaxisargs = Dict(:xrectzoom => false, :yrectzoom=>false, :titlevisible=>false, 
    :xgridvisible=>false, :ygridvisible=>false, 
    :xlabelvisible=>false, :ylabelvisible=>false, 
    :bottomspinevisible=>false, :leftspinevisible=>false, :topspinevisible=>false ,:rightspinevisible=>false,
    :xticklabelsvisible=>false, :yticklabelsvisible=>false, 
    :xticksvisible=>false, :yticksvisible=>false)
    
    fig = if isa(fig, Figure)
        my_layout = fig[1,1] = GridLayout()
    else
        my_layout = fig[] = GridLayout()
    end

    show_cbar = colorbar

    ax_im = ax_xz = ax_zy = nothing
    sl_z = nothing

    colorrange = Observable((0.0,1.0))

    if ndims(myim) > 2 && sz[3] > 1
        # grid_size = 3
        if ndims(myim) > 8
            error("Maximal number of dimensions exceeded.")
        end
        sl_x = Slider(fig[1,1,Bottom()], range = 1:sz[1], horizontal = true, startvalue = sz[1]/2+1)
        sl_y = Slider(fig[1,1,Right()], range = sz[2]:-1:1, horizontal = false, startvalue = sz[2]/2+1)
        sl_z = Slider(fig[1,2,Bottom()], range = 1:sz[3], horizontal = true, startvalue = sz[3]/2+1)

        sl_o = Tuple(Slider(fig[1:2,2+show_cbar+d], range = sz[3+d]:-1:1, horizontal = false, startvalue = 1) for d=1:ndims(myim)-3)
        pos_o = obs_from_sliders(sl_o)
        ax_im = Axis(fig[1,1], yreversed = true, alignmode = Inside(); hideaxisargs...)
        ax_im_xz = Axis(fig[2,1], yreversed = true, alignmode = Inside(); hideaxisargs...) # 
        ax_im_zy = Axis(fig[1,2], yreversed = true, alignmode = Inside(); hideaxisargs...) # 


        my_cmap = colormap_line(fig[3,1:3]; textsize=textsize) # 554 Mb
        if isnothing(pos_o)
            myim_xz = @lift(collect(get_slice(myim, (sl_x.value, $(sl_y.value), sl_z.value), (1,3))))
            myim_zy = @lift(collect(get_slice(myim, ( $(sl_x.value), sl_y.value, sl_z.value), (3,2))))
            myim_xy = @lift(collect(get_slice(myim, (sl_x.value, sl_y.value, $(sl_z.value)), (1,2))))
            position = @lift(get_pos(($(sl_x.value), $(sl_y.value), $(sl_z.value))))
        else
            myim_xz = @lift(collect(get_slice(myim, (sl_x.value, $(sl_y.value), sl_z.value, $(pos_o)...), (1,3))))
            myim_zy = @lift(collect(get_slice(myim, ( $(sl_x.value), sl_y.value, sl_z.value, $(pos_o)...), (3,2))))
            myim_xy = @lift(collect(get_slice(myim, (sl_x.value, sl_y.value, $(sl_z.value), $(pos_o)...), (1,2))))
            position = @lift(get_pos(($(sl_x.value), $(sl_y.value), $(sl_z.value), $(pos_o)...)))
        end
        im = heatmap!(fig[1,1], myim_xy, interpolate=false, highclip=:red, lowclip=:blue, colormap=my_cmap)
        xlims!(0,sz[1])
        ylims!(sz[2],0) # reverse y !
        crosshair(ax_im, sl_x, sl_y, (sz[1],sz[2]), color=color)

        im_xz = heatmap!(fig[2,1], myim_xz, interpolate=false, highclip=:red, lowclip=:blue, colormap=my_cmap)
        xlims!(0,sz[1])
        ylims!(sz[3],0) # no reverse
        crosshair(ax_im_xz, sl_x, sl_z, (sz[1],sz[3]), color=color)

        im_zy = heatmap!(fig[1,2], myim_zy, interpolate=false, highclip=:red, lowclip=:blue, colormap=my_cmap)
        xlims!(0,sz[3])
        ylims!(sz[1],0) # no reverse
        crosshair(ax_im_zy, sl_z, sl_y, (sz[3],sz[2]), color=color)

        colorrange[]=(
            min(to_value(im.attributes.colorrange)[1],to_value(im_zy.attributes.colorrange)[1],to_value(im_xz.attributes.colorrange)[1]),
            max.(to_value(im.attributes.colorrange)[2],to_value(im_zy.attributes.colorrange)[2],to_value(im_xz.attributes.colorrange)[2]))
    
        im.attributes.colorrange=colorrange
        im_xz.attributes.colorrange=colorrange
        im_zy.attributes.colorrange=colorrange

        ref_ax = [ax_im]

        txt = @lift(get_text(title, to_value($(position)), myim, aspects))
        my_label = Label(fig[2,2:2+show_cbar],txt, halign = :left, valign = :top, tellheight=false, textsize=textsize) # 

        if show_cbar
            # label = "Brightness", 
            cbar = Colorbar(fig[1,3], im, alignmode = Outside(), width=15, ticklabelspace=25f0)
            register_colorbar_scroll!(cbar, colorrange, to_value(colorrange))
        end
        r_ratio = Auto(aspects[3]*sz[3]/(aspects[2]*sz[2]))
        c_ratio = Auto(aspects[3]*sz[3]/(aspects[1]*sz[1]))

        colsize!(fig, 2, c_ratio)
        rowsize!(my_layout, 2, r_ratio)
        rowsize!(my_layout, 3, 30)

        register_panel_interactions!(ax_im_xz, sl_x, sl_z, sl_y, ref_ax)
        register_panel_interactions!(ax_im_zy, sl_z, sl_y, sl_x, ref_ax)
        register_panel_zoom_link!(ax_im, position, sz, ax_im_xz, ax_im_zy, aspects=aspects)
    else
        ax_im = Axis(fig[1,1], title=title, yreversed = true, alignmode = Inside(); hideaxisargs...)    
    
        sl_x = Slider(fig[1,1, Bottom()], range = 1:sz[1], horizontal = true, startvalue = sz[1]/2+1)
        sl_y = Slider(fig[1,1, Right()], range = sz[2]:-1:1, horizontal = false, startvalue = sz[2]/2+1)
        sl_o = Tuple(Slider(fig[1:2,1+d+show_cbar], range = sz[3+d]:-1:1, horizontal = false, startvalue = 1) for d=1:ndims(myim)-3)
        pos_o = obs_from_sliders(sl_o)
        # @time hidedecorations!(ax_im, grid = false);
        myim_xy, position = let 
            if isnothing(pos_o)
                myim_xy = collect(get_slice(myim, (sl_x.value, sl_y.value, 1), (1,2)))
                if ndims(myim) == 3
                    position = @lift(get_pos(($(sl_x.value), $(sl_y.value),1)))
                else
                    position = @lift(get_pos(($(sl_x.value), $(sl_y.value))))
                end
                myim_xy, position
            else
                myim_xy = @lift(collect(get_slice(myim, (sl_x.value, sl_y.value, 1, $(pos_o)...), (1,2)))) # is needed due to other positions
                position = @lift(get_pos(($(sl_x.value), $(sl_y.value), 1, $(pos_o)...)))
                myim_xy, position
            end
        end
        my_cmap = colormap_line(fig[3,1:2]; textsize=textsize)    
        im = heatmap!(myim_xy, interpolate=false, highclip=:red, lowclip=:blue, colormap=my_cmap)

        xlims!(0,sz[1])
        ylims!(sz[2],0) # reverse y !
        crosshair(ax_im, sl_x,sl_y, sz[1:2], color=color)

        colorrange[] = to_value(im.attributes.colorrange)    
        im.attributes.colorrange=colorrange

        ref_ax = [ax_im]
        txt = @lift(get_text(title, to_value($(position)), myim, aspects))
        my_label = Label(fig[2,1:1+show_cbar],txt, halign = :left, valign = :top, textsize=textsize)

        if show_cbar
            # label = "Brightness", 
            cbar = Colorbar(fig[1,2], im, alignmode = Outside(), width=20, ticklabelspace=30f0)
            register_colorbar_scroll!(cbar, colorrange, to_value(colorrange))
        end

        register_panel_zoom_link!(ax_im, position, sz, nothing, nothing, aspects=aspects)

    end

    register_panel_interactions!(ax_im, sl_x, sl_y, sl_z, ref_ax)

    abs_zoom = get_max_zoom(sz.*aspects, ax_im, ax_zy, ax_xz)
    zoom_all(abs_zoom, nothing, ax_im, ax_xz, ax_zy, aspects, sz .÷ 2 .+1, sz)

    return fig # , grid_size
end

function ortho!(myim; preferred_size = 600, title = "Image", color=:red, markersize = 40.0, aspects=ones(ndims(myim)), colorbar=true, textsize=24)
    sz = size(myim) .* aspects
    if length(sz) > 2 && sz[3] != 1
        r_ratio = sz[3]/sz[2] # these ratios will be used in the ortho! function to specify the relative size
        c_ratio = sz[3]/sz[1]
        sz = sz[1:3]
        fak = preferred_size ./ max(sz...)
        sz3 = max(sz[3],160/fak)
        res = fak .* (sz[1]+sz3, sz[2]+sz[3]) ./ 2 .+ (80*colorbar,70)
    else
        fak = preferred_size ./ max(sz...)
        res = fak .* (sz[1], sz[2])  .+ (30*colorbar,160) # additional space for cbar and text
    end
    fig = Figure(resolution=res)
    ortho!(fig, myim; title = title,  color=color, markersize = markersize, aspects=aspects, colorbar=colorbar, textsize=textsize)
    return fig
end

function ortho(myim; kwargs...)
    myfig = ortho!(myim; kwargs...)
    ns = GLMakie.Screen()
    GLMakie.display(ns, myfig)
end

function test_ortho_view(;aspects=(1,1,1,1))
    set_theme!(theme_black())
    obj = rand(10,20,30,40) # testimage("simple_3d_ball.tif")
    ortho!(obj, aspects=aspects)
end

end
