
import modules as m

def test():
	
	svg_doc = m.svg.Svg(dpi=72)

	# text_box_height = svg_doc.text(
	# 	text = 'AcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumannii',
	# 	x = 39.5,
	# 	y = 659,
	# 	text_box_width = (172.11 - 39.5) - (179.46 - 173.11),
	# 	font_size = 12,
	# 	font_family = 'Barlow',
	# 	font_ttf = 'Barlow-Regular.ttf'
	# )
	# print(f'text_box_height: {text_box_height}')

	text_box_height = svg_doc.text(
		text = 'AcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumanniiAcinetobacterbaumannii',
		x = 39.5,
		y = 659,
		text_box_width = (172.11 - 39.5) - (179.46 - 173.11),
		font_size = 12,
		font_family = 'Barlow',
		font_ttf = 'Barlow-Regular.ttf',
		bullet = '•'
	)

	svg_doc.export(filetype='svg')
	svg_doc.export(filetype='png')