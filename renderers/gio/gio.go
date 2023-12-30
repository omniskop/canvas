package gio

import (
	"image"
	"image/color"
	"math"

	"gioui.org/f32"
	"gioui.org/layout"
	"gioui.org/op"
	"gioui.org/op/clip"
	"gioui.org/op/paint"
	"github.com/tdewolff/canvas"
)

type Gio struct {
	ops            *op.Ops
	width, height  float64
	xScale, yScale float64
	dimensions     layout.Dimensions
}

// New returns a Gio renderer of fixed size.
func New(gtx layout.Context, width, height float64) *Gio {
	dimensions := layout.Dimensions{Size: image.Point{int(width + 0.5), int(height + 0.5)}}
	return &Gio{
		ops:        gtx.Ops,
		width:      width,
		height:     height,
		xScale:     1.0,
		yScale:     1.0,
		dimensions: dimensions,
	}
}

// NewContain returns a Gio renderer that fills the constraints either horizontally or vertically, whichever is met first.
func NewContain(gtx layout.Context, width, height float64) *Gio {
	xScale := float64(gtx.Constraints.Max.X-gtx.Constraints.Min.X) / width
	yScale := float64(gtx.Constraints.Max.Y-gtx.Constraints.Min.Y) / height
	if yScale < xScale {
		xScale = yScale
	} else {
		yScale = xScale
	}

	dimensions := layout.Dimensions{Size: image.Point{int(width*xScale + 0.5), int(height*yScale + 0.5)}}
	return &Gio{
		ops:        gtx.Ops,
		width:      width,
		height:     height,
		xScale:     xScale,
		yScale:     yScale,
		dimensions: dimensions,
	}
}

// NewStretch returns a Gio renderer that stretches the view to fit the constraints.
func NewStretch(gtx layout.Context, width, height float64) *Gio {
	xScale := float64(gtx.Constraints.Max.X-gtx.Constraints.Min.X) / width
	yScale := float64(gtx.Constraints.Max.Y-gtx.Constraints.Min.Y) / height

	dimensions := layout.Dimensions{Size: image.Point{int(width*xScale + 0.5), int(height*yScale + 0.5)}}
	return &Gio{
		ops:        gtx.Ops,
		width:      width,
		height:     height,
		xScale:     xScale,
		yScale:     yScale,
		dimensions: dimensions,
	}
}

// Dimensions returns the dimensions for Gio.
func (r *Gio) Dimensions() layout.Dimensions {
	return r.dimensions
}

// Size returns the size of the canvas in millimeters.
func (r *Gio) Size() (float64, float64) {
	return r.width, r.height
}

func (r *Gio) point(p canvas.Point) f32.Point {
	return f32.Point{float32(r.xScale * p.X), float32(r.yScale * (r.height - p.Y))}
}

func (r *Gio) renderPath(path *canvas.Path, fill canvas.Paint) {
	path = path.ReplaceArcs()

	p := clip.Path{}
	p.Begin(r.ops)
	for scanner := path.Scanner(); scanner.Scan(); {
		switch scanner.Cmd() {
		case canvas.MoveToCmd:
			p.MoveTo(r.point(scanner.End()))
		case canvas.LineToCmd:
			p.LineTo(r.point(scanner.End()))
		case canvas.QuadToCmd:
			p.QuadTo(r.point(scanner.CP1()), r.point(scanner.End()))
		case canvas.CubeToCmd:
			p.CubeTo(r.point(scanner.CP1()), r.point(scanner.CP2()), r.point(scanner.End()))
		case canvas.ArcToCmd:
			// TODO: ArcTo
			p.LineTo(r.point(scanner.End()))
		case canvas.CloseCmd:
			p.Close()
		}
	}

	shape := clip.Outline{p.End()}
	defer shape.Op().Push(r.ops).Pop()

	if fill.IsColor() {
		paint.Fill(r.ops, toNRGBA(fill.Color))
	} else if fill.IsGradient() {
		if g, ok := fill.Gradient.(*canvas.LinearGradient); ok && len(g.Stops) == 2 {
			linearGradient := paint.LinearGradientOp{}
			linearGradient.Stop1 = r.point(g.Start)
			linearGradient.Stop2 = r.point(g.End)
			linearGradient.Color1 = toNRGBA(g.Stops[0].Color)
			linearGradient.Color2 = toNRGBA(g.Stops[1].Color)
			linearGradient.Add(r.ops)
			paint.PaintOp{}.Add(r.ops)
		}
	}
}

/*
rx, ry, rot, large, sweep := scanner.Arc()
			rot = (rot / 360) * 2 * math.Pi // convert from degrees to radiants
			center, _, radiusDelta := convertEllipse(scanner.Start(), scanner.End(), rx, ry, rot, large, sweep)

			centerToFocal := math.Sqrt(rx*rx - ry*ry)
			f1 := movePointInDirection(center, rot, centerToFocal)
			f2 := movePointInDirection(center, rot, -centerToFocal)

			p.ArcTo(r.point(f2), r.point(f1), float32(-radiusDelta))
*/

// convertEllipse takes parameters describing an ellipse in 'endpoint parameterization' and converts it to a 'center parameterization'.
// This was mainly taken from https://www.w3.org/TR/SVG/implnote.html#ArcConversionEndpointToCenter and https://gist.github.com/balint42/fdb1d7d2e16fe11ac785.
// The rotation must be in *radiants*.
func convertEllipse(start canvas.Point, end canvas.Point, rx, ry, rot float64, large, sweep bool) (center canvas.Point, startRadius, deltaRadius float64) {
	// source: https://www.w3.org/TR/SVG/implnote.html#ArcConversionEndpointToCenter

	// function for calculating angle between two vectors
	angle := func(u, v canvas.Point) float64 {
		var sign float64 = -1
		if (u.X*v.Y - u.Y*v.X) > 0 {
			sign = 1
		}
		return sign * math.Acos(
			(u.X*v.X+u.Y*v.Y)/
				(math.Sqrt(u.X*u.X+u.Y*u.Y)*math.Sqrt(u.X*u.X+u.Y*u.Y)),
		)
	}

	// sanitize input
	rot = math.Mod(rot, math.Pi*2)
	rx = math.Abs(rx)
	ry = math.Abs(ry)

	// Step 1:
	// Calculate vector from end point to middle between start end end.
	middle := canvas.Point{(start.X - end.X) / 2, (start.Y - end.Y) / 2}
	// Rotate that vector by -rot.
	cosRot := math.Cos(rot) // x component of rotation
	sinRot := math.Sin(rot) // y component of rotation
	x := cosRot*middle.X + sinRot*middle.Y
	y := -1*sinRot*middle.X + cosRot*middle.Y
	// Step 2: calculate center point relative to the middle point.
	var rx2 = rx * rx
	var ry2 = ry * ry
	var x2 = x * x
	var y2 = y * y
	var sign float64 = 1
	if large == sweep {
		sign = -1
	}
	var fr = sign * math.Sqrt(
		(rx2*(ry2-y2)-ry2*x2)/
			(rx2*y2+ry2*x2),
	)
	var xt = fr * rx * y / ry
	var yt = -1 * fr * ry * x / rx

	// Step 3: Reverse rotation and convert relative to absolute coordinates.
	var cx = cosRot*xt - sinRot*yt + (start.X+end.X)/2
	var cy = sinRot*xt + cosRot*yt + (start.Y+end.Y)/2

	// Step 4: calculate angles
	var vt = canvas.Point{X: (x - xt) / rx, Y: (y - yt) / ry}
	var phi1 = angle(canvas.Point{X: 1, Y: 0}, vt)
	var phiD = math.Mod(angle(vt, canvas.Point{X: (-x - xt) / rx, Y: (-y - yt) / ry}), math.Pi*2)

	return canvas.Point{X: cx, Y: cy}, phi1, phiD
}

func movePointInDirection(p canvas.Point, angle float64, distance float64) canvas.Point {
	return canvas.Point{
		p.X + math.Cos(angle)*distance,
		p.Y * math.Sin(angle) * distance,
	}
}

// RenderPath renders a path to the canvas using a style and a transformation matrix.
func (r *Gio) RenderPath(path *canvas.Path, style canvas.Style, m canvas.Matrix) {
	if style.HasFill() {
		r.renderPath(path.Transform(m), style.Fill)
	}

	if style.HasStroke() {
		if style.IsDashed() {
			path = path.Dash(style.DashOffset, style.Dashes...)
		}
		path = path.Stroke(style.StrokeWidth, style.StrokeCapper, style.StrokeJoiner, canvas.Tolerance)
		r.renderPath(path.Transform(m), style.Stroke)
	}
}

// RenderText renders a text object to the canvas using a transformation matrix.
func (r *Gio) RenderText(text *canvas.Text, m canvas.Matrix) {
	text.RenderAsPath(r, m, 0.0)
}

// RenderImage renders an image to the canvas using a transformation matrix.
func (r *Gio) RenderImage(img image.Image, m canvas.Matrix) {
	paint.NewImageOp(img).Add(r.ops)
	m = canvas.Identity.Scale(r.xScale, r.yScale).Mul(m)
	m = m.Translate(0.0, float64(img.Bounds().Max.Y))
	trans := op.Affine(f32.NewAffine2D(
		float32(m[0][0]), -float32(m[0][1]), float32(m[0][2]),
		-float32(m[1][0]), float32(m[1][1]), float32(r.yScale*r.height-m[1][2]),
	)).Push(r.ops)
	paint.PaintOp{}.Add(r.ops)
	trans.Pop()
}

func toNRGBA(col color.Color) color.NRGBA {
	r, g, b, a := col.RGBA()
	if a == 0 {
		return color.NRGBA{}
	}
	r = (r * 0xffff) / a
	g = (g * 0xffff) / a
	b = (b * 0xffff) / a
	return color.NRGBA{R: uint8(r >> 8), G: uint8(g >> 8), B: uint8(b >> 8), A: uint8(a >> 8)}
}
