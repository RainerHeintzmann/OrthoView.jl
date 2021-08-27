# cd("C:\\Users\\pi96doc\\Documents\\Programming\\Julia\\OrthoView.jl\\")
# ]activate .
## interactive example:
module OrthoView
export ortho!, test_ortho_view

using GLMakie
using GLMakie.FileIO
using Printf
using ColorTypes  # for RGB type
using ColorSchemes

function get_crosshair_xs(px, sx, gap)
    xs = [0f0,px-gap,px+gap,sx,px,px,px,px]
    return xs
end
function get_crosshair_ys(py, sy, gap)
    gap = Float32(to_value(gap))
    ys = [py,py,py,py,0f0,py-gap,py+gap,sy]
    return ys
end

function crosshair(sl_x,sl_y, sz; markersize=0.7, color=:red)
    if isa(markersize, Slider)
        cross_xs = @lift(get_crosshair_xs($(sl_x.value)-0.5f0, sz[1], Float32($(markersize.value))))
        cross_ys = @lift(get_crosshair_ys($(sl_y.value)-0.5f0, sz[2], Float32($(markersize.value))))
    else
        gap = Float32(to_value(markersize))
        cross_xs = @lift(get_crosshair_xs($(sl_x.value)-0.5f0, sz[1], gap))
        cross_ys = @lift(get_crosshair_ys($(sl_y.value)-0.5f0, sz[2], gap))
    end
    linesegments!(cross_xs, cross_ys, color=color, linewidth=3)
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
    lx,ly = to_value(ax_xy.limits)
    sx,sy = ax_xy.layoutobservables.computedbbox.val.widths  # The size in pixels that is available
    if isnothing(lx)
        return Consume(true)
    end
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

function zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
    if !isnothing(event)
        zoom_dir = to_value(event.y)
        zoom_max = get_max_zoom(sz .* aspects, ax, ax_zy, ax_xz)
        abs_zoom = min(abs_zoom * 1.1 ^ -zoom_dir, zoom_max)
    end
    do_zoom(ax, ctr, sz, abs_zoom, aspects, (1,-2))
    do_zoom(ax_xz, ctr, sz, abs_zoom, aspects, (1,-3))
    do_zoom(ax_zy, ctr, sz, abs_zoom, aspects, (3,-2))
    return abs_zoom
end

function register_panel_zoom_link!(ax, ctr, sz, ax_xz=nothing, ax_zy=nothing; aspects=ones(length(sz)))

    abs_zoom = get_max_zoom(sz.*aspects, ax, ax_zy, ax_xz)
    zoom_all(abs_zoom, nothing, ax, ax_xz, ax_zy, aspects, ctr, sz)
    register_interaction!(ax, :my_scroll_interaction) do event::ScrollEvent, axis
        abs_zoom = zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
        return Consume(true) # for now prevent the zoom
    end
    if !isnothing(ax_xz)
        register_interaction!(ax_xz, :my_scroll_interaction) do event::ScrollEvent, axis
            abs_zoom = zoom_all(abs_zoom, event, ax, ax_xz, ax_zy, aspects, ctr, sz)
            return Consume(true) # for now prevent the zoom
        end
    end
    if !isnothing(ax_zy)
        register_interaction!(ax_zy, :my_scroll_interaction) do event::ScrollEvent, axis
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

function register_panel_interactions!(ax, sl_x, sl_y, sl_z, ref_ax; key_buffer="")
    deregister_interaction!(ax, :rectanglezoom)
    sz = (max(to_value(sl_x.range)[1],to_value(sl_x.range)[end]), max(to_value(sl_y.range)[1],to_value(sl_y.range)[end])) # (to_value(ax.limits)[1][2] , to_value(ax.limits)[2][2])
    ctr = sz .÷ 2 .+ 1
    # pos = Node([Point2f0(0)])

    register_interaction!(ax, :my_mouse_interaction) do event::MouseEvent, axis
        if event.type === MouseEventTypes.leftclick || event.type === MouseEventTypes.leftdragstart || 
            event.type === MouseEventTypes.leftdragstop || event.type === MouseEventTypes.leftdrag || 
            event.type === MouseEventTypes.middledragstart || event.type === MouseEventTypes.middledragstop || event.type === MouseEventTypes.middledrag
            set_close_to!(sl_x, event.data[1] + 0.5f0)
            set_close_to!(sl_y, event.data[2] + 0.5f0)
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

function get_slice(data, pos, dims=(1,2), min_max_val=extrema(data)) 
    idx = Tuple((d in dims) ? Colon() : pos[d] for d = 1:length(pos))
    if dims[1] > dims[2]
        res = transpose(@view data[idx...])
    else
        res = @view data[idx...]
    end
    if isa(res[1],Complex)
        res = abs2.(res)
    end
    # res = map(apply_min_max_gamma(min_max_val[1], min_max_val[2], 1.0f0), res)
    return res
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

function ortho!(fig, myim; title = "Image", color=:red, markersize = 40.0, aspects=ones(ndims(myim)), colorbar=true)
    sz = size(myim)

    fig = if isa(fig, Figure)
        fig[1,1] = GridLayout()
    else
        fig[] = GridLayout()
    end

    show_cbar = colorbar

    ax_im = ax_xz = ax_zy = nothing
    sl_z = nothing

    ColOptions = ["greys","greens","reds", "blues","RdBu","heat", "viridis", "plasma", "magma", "inferno"]

    if ndims(myim) > 2 && sz[3] > 1
        # grid_size = 3
        sl_x = Slider(fig[1,1,Bottom()], range = 1:sz[1], horizontal = true, startvalue = sz[1]/2+1)
        sl_y = Slider(fig[1,1,Right()], range = sz[2]:-1:1, horizontal = false, startvalue = sz[2]/2+1)
        sl_z = Slider(fig[1,2,Bottom()], range = 1:sz[3], horizontal = true, startvalue = sz[3]/2+1)
        if ndims(myim) > 8
            error("Maximal number of dimensions exceeded.")
        end

        sl_o = Tuple(Slider(fig[1:2,2+show_cbar+d], range = sz[3+d]:-1:1, horizontal = false, startvalue = 1) for d=1:ndims(myim)-3)
        pos_o = obs_from_sliders(sl_o)
        ax_im_xz = Axis(fig[2,1], yreversed = true, alignmode = Inside(), xrectzoom=false, yrectzoom=false) # 
        hidedecorations!(ax_im_xz, grid = false); ax_im_xz.xrectzoom = false; ax_im_xz.yrectzoom = false
        # myim_xz = @lift(myim[:,round(Int,$(sl_y.value)),:])
        if isnothing(pos_o)
            myim_xz = @lift(get_slice(myim, (sl_x.value, $(sl_y.value), sl_z.value), (1,3)))
            myim_zy = @lift(get_slice(myim, ( $(sl_x.value), sl_y.value, sl_z.value), (3,2)))
            myim_xy = @lift(get_slice(myim, (sl_x.value, sl_y.value, $(sl_z.value)), (1,2)))
            position = @lift(get_pos(($(sl_x.value), $(sl_y.value), $(sl_z.value))))
        else
            myim_xz = @lift(get_slice(myim, (sl_x.value, $(sl_y.value), sl_z.value, $(pos_o)...), (1,3)))
            myim_zy = @lift(get_slice(myim, ( $(sl_x.value), sl_y.value, sl_z.value, $(pos_o)...), (3,2)))
            myim_xy = @lift(get_slice(myim, (sl_x.value, sl_y.value, $(sl_z.value), $(pos_o)...), (1,2)))
            position = @lift(get_pos(($(sl_x.value), $(sl_y.value), $(sl_z.value), $(pos_o)...)))
        end
        im_xz = image!(myim_xz, interpolate=false)
        xlims!(0,sz[1])
        ylims!(sz[3],0) # no reverse
        crosshair(sl_x,sl_z, (sz[1],sz[3]), color=color)

        ax_im_zy = Axis(fig[1,2], yreversed = true, alignmode = Inside(), xrectzoom=false, yrectzoom=false) # 
        hidedecorations!(ax_im_zy, grid = false); ax_im_zy.xrectzoom = false; ax_im_zy.yrectzoom = false
        im_zy = image!(myim_zy, interpolate=false)
        xlims!(0,sz[3])
        ylims!(sz[1],0) # no reverse
        crosshair(sl_z,sl_y, (sz[3],sz[2]), color=color)

        #  aspect = DataAspect(), 
        ax_im = Axis(fig[1,1], yreversed = true, alignmode = Inside(), xrectzoom=false, yrectzoom=false)
        hidedecorations!(ax_im, grid = false); ax_im.xrectzoom = false; ax_im.yrectzoom = false
        im = image!(myim_xy, interpolate=false)
        #xlims!(4,5)
        #ylims!(4,5)
        # ax_im.aspect = DataAspect()
        xlims!(0,sz[1])
        ylims!(sz[2],0) # reverse y !
        crosshair(sl_x,sl_y, (sz[1],sz[2]), color=color)

        # linkxaxes!(ax_im, ax_im_xz) # to ensure the scrolling behaves correctly
        # linkyaxes!(ax_im, ax_im_zy)

        ref_ax = [ax_im]

        #txt = @lift(get_text(title, to_value($(position)), myim, aspects))
        #ax_txt = Textbox(fig[2,2:2+show_cbar], placeholder = txt, width = 300)
    
        ax_txt = Axis(fig[2,2:2+show_cbar], backgroundcolor=:white, xzoomlock=true, yzoomlock=true, xpanlock=true, ypanlock=true, xrectzoom=false, yrectzoom=false) # title=title, 
        hidedecorations!(ax_txt, grid = false); ax_im.xrectzoom = false; ax_im.yrectzoom = false

        r_ratio = Auto(aspects[3]*sz[3]/(aspects[2]*sz[2]))
        c_ratio = Auto(aspects[3]*sz[3]/(aspects[2]*sz[2]))
        rowsize!(fig, 2, r_ratio)
        colsize!(fig, 2, c_ratio)

        # @show ax_txt.layoutobservables.computedbbox.val.widths
        # if ax_txt.layoutobservables.computedbbox.val.widths[1] < 350 # strange why this does not correspond to the number to set to...
        #     colsize!(fig, 2, 180)
        # end
        # if ax_txt.layoutobservables.computedbbox.val.widths[2] < 130
        #     rowsize!(fig, 2, 130)
        # end
        # @show ax_txt.layoutobservables.computedbbox.val.widths
        if show_cbar
            cbar = Colorbar(fig[1,3], im, label = "Brightness", alignmode = Outside())
            cbar.width = 20
        end

        register_panel_interactions!(ax_im_xz, sl_x, sl_z, sl_y, ref_ax)
        register_panel_interactions!(ax_im_zy, sl_z, sl_y, sl_x, ref_ax)
        register_panel_zoom_link!(ax_im, position, sz, ax_im_xz, ax_im_zy, aspects=aspects)
        menu = Menu(fig[3,1], options = ColOptions)
        on(menu.selection) do s
            im.colormap = s
            im_xz.colormap = s
            im_zy.colormap = s
         end
    else
        ax_im = Axis(fig[1,1], title=title, yreversed = true, alignmode = Inside(), xrectzoom=false, yrectzoom=false)
        sl_x = Slider(fig[1,1, Bottom()], range = 1:sz[1], horizontal = true, startvalue = sz[1]/2+1)
        sl_y = Slider(fig[1,1, Right()], range = sz[2]:-1:1, horizontal = false, startvalue = sz[2]/2+1)
        sl_o = Tuple(Slider(fig[1:2,1+d+show_cbar], range = sz[3+d]:-1:1, horizontal = false, startvalue = 1) for d=1:ndims(myim)-3)
        pos_o = obs_from_sliders(sl_o)
        hidedecorations!(ax_im, grid = false); ax_im.xrectzoom = false; ax_im.yrectzoom = false
        if isnothing(pos_o)
            myim_xy = get_slice(myim, (sl_x.value, sl_y.value, 1), (1,2))
            if ndims(myim) == 3
                position = @lift(get_pos(($(sl_x.value), $(sl_y.value),1)))
            else
                position = @lift(get_pos(($(sl_x.value), $(sl_y.value))))
            end
        else
            myim_xy = @lift(get_slice(myim, (sl_x.value, sl_y.value, 1, $(pos_o)...), (1,2)))
            position = @lift(get_pos(($(sl_x.value), $(sl_y.value), 1, $(pos_o)...)))
        end
        im = image!(myim_xy, interpolate=false)
        xlims!(0,sz[1])
        ylims!(sz[2],0) # reverse y !
        crosshair(sl_x,sl_y, sz[1:2],color=color)
        ref_ax = [ax_im]
        #txt = @lift(get_text(title, to_value($(position)), myim, aspects))
        #ax_txt = Textbox(fig[2,1:1+show_cbar], placeholder = txt, width = 300)

        ax_txt = Axis(fig[2,1:1+show_cbar], backgroundcolor=:white, xzoomlock=true, yzoomlock=true, xpanlock=true, ypanlock=true, xrectzoom=false, yrectzoom=false) # title=title, 
        rowsize!(fig, 2, 180)
        hidedecorations!(ax_txt, grid = false); ax_im.xrectzoom = false; ax_im.yrectzoom = false

        if show_cbar
            cbar = Colorbar(fig[1,2], im, label = "Brightness", alignmode = Outside())
            cbar.width = 20
        end

        register_panel_zoom_link!(ax_im, position, sz, nothing, nothing, aspects=aspects)
        menu = Menu(fig[3,1], options = ColOptions)
        on(menu.selection) do s
            im.colormap = s
         end     
    end

    register_panel_interactions!(ax_im, sl_x, sl_y, sl_z, ref_ax)
    xlims!(ax_txt,-4,100); ylims!(ax_txt,-4,100)
    txt = @lift(get_text(title, to_value($(position)), myim, aspects))
    text!(ax_txt, txt, position = Point(0,99), color = :black, align = (:left, :top), justification = :left)

    abs_zoom = get_max_zoom(sz.*aspects, ax_im, ax_zy, ax_xz)
    zoom_all(abs_zoom, nothing, ax_im, ax_xz, ax_zy, aspects, sz .÷ 2 .+1, sz)

    return fig # , grid_size
end

function ortho!(myim; preferred_size = 600, title = "Image", color=:red, markersize = 40.0, aspects=ones(ndims(myim)), colorbar=true)
    sz = size(myim) .* aspects
    fak = preferred_size ./ max(sz...)
    if length(sz) > 2 && sz[3] != 1
        sz3 = max(fak .* sz[3],210)
        res = fak .* (sz[1]+sz3, sz[2]+sz[3])
    else
        res = fak .* (sz[1], sz[2])  .+ (30*colorbar,220) # additional space for cbar and text
    end
    fig = Figure(resolution=res)
    ortho!(fig, myim; title = title,  color=color, markersize = markersize, aspects=aspects, colorbar=colorbar)
    return fig
end

function test_ortho_view(;aspects=(1,1,1,1))
    set_theme!(theme_black())
    obj = rand(10,20,30,40) # testimage("simple_3d_ball.tif")
    ortho!(obj, aspects=aspects)
end

end
