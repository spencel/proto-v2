
import drawsvg
from PIL import (
	Image as PilImage,
	ImageDraw as PilImageDraw,
	ImageFont as PilImageFont
)


class Svg:

	document_sizes = {
		'letter_size': {
			'width_in': 8.5, # inches
			'height_in': 11 # inches
		}
	}

	def __init__(self,
		width:int = 0,  # px
		height:int = 0, # px
		dpi:int = 72,   # px
			# [pixel qty] = DPI / 72
			# 72 DPI: 1 point = 1 pixel
			# 144 DPI: 1 point = 2 pixels
			# 288 DPI: 1 point = 4 pixels
		document_size:str = 'letter_size'
	):
		if width or height:
			self.width = width
			self.height = height
		elif dpi:
			self.width = dpi * self.document_sizes[document_size]['width_in']
			self.height = dpi * self.document_sizes[document_size]['height_in']
		else:
			raise Exception('Error: Svg class incorrectly instantiated.')
		self.dpi = dpi
		self.document_size = document_size

		self.drawing = drawsvg.Drawing(
			width = self.width,
			height = self.height,
			origin = (0, 0),
			context = None, # drawsvg.types.Context
			animation_config = None,
			id_prefix= 'd',
		)


	def append_raw_svg(self, raw_svg):
		self.drawing.append(drawsvg.Raw(raw_svg))


	def text(self,
		text:str = 'Text Element',
		x:int = 0,
		y:int = 0,
		text_box_width = None,
		font_size:int = 12, # pt, 8.63 px
		font_family:str = 'Arial',
		font_ttf:str = '',
		bullet:str = '',
		center:bool = False,
		leading:int = 1.2, # em, 5.77 px
		line_offset:int = 0,
		path = None,
		start_offset = None,
		path_args = None,
		tspan_args = None,
		is_cairo_fix:bool = True
	):
		if text_box_width and not font_ttf:
			raise Exception('Error: Font ttf was not provided for text box.')
		
		if bullet:
			self.drawing.append(drawsvg.Text(
				bullet,
				font_size,
				x,
				y,
				font_family = font_family
			))
			# Get width of bullet
			dum_image = PilImage.new('RGB', (1, 1))
			dum_draw = PilImageDraw.Draw(dum_image)
			font = PilImageFont.truetype(font_ttf, font_size)
			textbbox = dum_draw.textbbox((0, 0), bullet, font=font)
			dum_width = textbbox[2] - textbbox[0]
			# Adjust new textbox width
			x += dum_width + 6.88
			text_box_width -= dum_width + 6.88
		
		dum_image = PilImage.new('RGB', (1, 1))
		dum_draw = PilImageDraw.Draw(dum_image)
		font = PilImageFont.truetype(font_ttf, font_size)
		text_lines = []
		total_height = 0
		if not text_box_width:
			text_lines = [text]
		else:
			temp_text = ''
			for i in range(len(text)):
				temp_text += text[i]
				textbbox = dum_draw.textbbox((0, 0), temp_text, font=font)
				dum_width = textbbox[2] - textbbox[0]
				if dum_width > text_box_width:
					# Add text line to list without this char
					text_lines.append(temp_text[:-1])
					# Calculate new textbox height
					dum_height = (textbbox[3] - textbbox[1]) * (8.63 / 8)
					print(f'dum_height: {dum_height}')
					total_height += dum_height + leading * (5.77 / 1.2)
					# Start new text line with this char
					temp_text = temp_text[-1]
			if len(temp_text) > 1:
				text_lines.append(temp_text)
			total_height -= leading * (5.77 / 1.2)
				
		self.drawing.append(drawsvg.Text(
			text_lines,
			font_size,
			x,
			y,
			font_family = font_family,
			center = center,
			line_height = leading,
			line_offset = line_offset,
			path = path,
			start_offset = start_offset,
			path_args = path_args,
			tspan_args = tspan_args,
			cairo_fix = is_cairo_fix
		))

		return total_height
	

	def export(self,
		fpath:str = '',
		filetype:str = 'svg'
	):
		match filetype:
			case 'svg':
				if not fpath:
					fpath = 'output.svg'
				self.drawing.save_svg(fpath)
			case 'png':
				if not fpath:
					fpath = 'output.png'
				self.drawing.save_png(fpath)
