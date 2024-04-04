
import re

class DataNormalizer():


    @staticmethod
    def normalize_date(raw_date, pattern="YYYY-MM-DD") -> str:

        year = "0000"
        month = "00"
        day = "00"

        if raw_date == "":
            return f"{year}-{month}-{day}"

        # Normalize character case
        raw_date = raw_date.lower()

        MONTH_NAME_MAP = {
            "january": "01",
            "jan": "01",
            "february": "02",
            "feb": "02",
            "march": "03",
            "mar": "03",
            "april": "04",
            "apr": "04",
            "may": "05",
            "june": "06",
            "jun": "06",
            "july": "07",
            "jul": "07",
            "august": "08",
            "aug": "08",
            "september": "09",
            "sep": "09",
            "october": "10",
            "oct": "10",
            "november": "11",
            "nov": "11",
            "december": "12",
            "dec": "12"
        }
        MONTHS_REGEX = "|".join(MONTH_NAME_MAP.keys())

        def check(month, day, raw_date):
            month = int(month)
            day = int(day)
            # Swap
            if (
                # Month is day
                month > 12
                and month <= 31
                # Day is month
                and day <= 12
            ):
                new_day = month
                month = day
                day = new_day
            # Bad quantities
            elif (
                month > 31
                or day > 31
                or (
                    month > 12
                    and day > 12
                )
            ):
                raise Exception(f"Unexpected date format: {raw_date}")

            if month < 10:
                month = "0" + str(month)
            if day < 10:
                day = "0" + str(day)

            return month, day


        # Only year provided #######################

        # YYYY
        if re.search("^\d{4}$", raw_date):
            year = raw_date

        # YY
        # Assume it's a year from the 1900s
        if re.search("^\d{2}$", raw_date):
            year = "19" + raw_date
            raise Exception(f"Unexpected date format: {raw_date}")

        # Only year and month provided ################

        # YYYY/MM
        if re.search("^\d{4}(/|-)\d{2}$", raw_date):
            year = raw_date[:4]
            month = raw_date[5:]

        # MM/YYYY
        if re.search("^\d{2}(/|-)\d{4}$", raw_date):
            year = raw_date[:2]
            month = raw_date[3:]

        # YYYY/M
        if re.search("^\d{4}(/|-)\d{1}$", raw_date):
            year = raw_date[:4]
            month = "0" + raw_date[5:]

        # M/YYYY
        if re.search("^\d{1}(/|-)\d{4}$", raw_date):
            year = raw_date[2:]
            month = "0" + raw_date[:1]

        # YY/MM
        if re.search("^\d{2}(/|-)\d{2}$", raw_date):
            year = "19" + raw_date[:2]
            month = raw_date[3:]
            raise Exception(f"Unexpected date format: {raw_date}")

        # MM/YY
        if re.search("^\d{1}(/|-)\d{4}$", raw_date):
            year = "19" + raw_date[3:]
            month = raw_date[:2]
            raise Exception(f"Unexpected date format: {raw_date}")

        # YY/M
        if re.search("^\d{2}(/|-)\d{1}$", raw_date):
            year = "19" + raw_date[:2]
            month = "0" + raw_date[3:]
            raise Exception(f"Unexpected date format: {raw_date}")

        # M/YY
        if re.search("^\d{1}(/|-)\d{2}$", raw_date):
            year = "19" + raw_date[2:]
            month = "0" + raw_date[:1]
            raise Exception(f"Unexpected date format: {raw_date}")

        # Month is named ##############################

        # M.../YYYY
        search_results = re.search(
            "^(" + MONTHS_REGEX + ")(/|-)(\d{4})$", raw_date
        )
        if search_results:
            month_name = search_results.group(1)
            month = MONTH_NAME_MAP[month_name.lower()]
            year = search_results.group(3)

        # M.../YY
        search_results = re.search(
            "^(" + MONTHS_REGEX + ")(/|-)(\d{2})$", raw_date
        )
        if search_results:
            month_name = search_results.group(1)
            month = MONTH_NAME_MAP[month_name.lower()]
            year = search_results.group(3)
            raise Exception(f"Unexpected date format: {raw_date}")

        # DD/M.../YYYY
        search_results = re.search(
            "^(\d{2})(/|-)(" + MONTHS_REGEX + ")(/|-)(\d{4})$", raw_date
        )
        if search_results:
            day = search_results.group(1)
            month_name = search_results.group(3)
            month = MONTH_NAME_MAP[month_name.lower()]
            year = search_results.group(5)

        # D/M.../YYYY
        search_results = re.search(
            "^(\d{1})(/|-)(" + MONTHS_REGEX + ")(/|-)(\d{2})$", raw_date
        )
        if search_results:
            day = "0" + search_results.group(1)
            month_name = search_results.group(3)
            month = MONTH_NAME_MAP[month_name.lower()]
            year = search_results.group(5)


        # Year, month, and day provided ################

        # YYYY/MM/DD
        if re.search("^\d{4}(/|-)\d{2}(/|-)\d{2}$", raw_date):
            year = raw_date[:4]
            month = int(raw_date[5:7])
            day = int(raw_date[8:])
            month, day = check(month, day, raw_date)

        # ??/??/YYYY
        # MM/DD/YYYY
        # DD/MM/YYYY
        if re.search("^\d{2}(/|-)\d{2}(/|-)\d{4}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = raw_date[6:]
            month, day = check(month, day, raw_date)

        # ??/??/??
        # YY/MM/DD
        # MM/DD/YY
        # DD/MM/YY
        if re.search("^\d{2}(/|-)\d{2}(/|-)\d{2}}$", raw_date):
            month = int(raw_date[:2])
            day = int(raw_date[3:5])
            year = "19" + raw_date[6:8]
            month, day = check(month, day, raw_date)
            raise Exception(f"Unexpected date format: {raw_date}")

        # ?/??/YYYY
        # M/DD/YYYY
        # D/MM/YYYY
        if re.search("^\d{1}(/|-)\d{2}(/|-)\d{4}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = raw_date[6:]
            month, day = check(month, day, raw_date)

        # ??/?/YYYY
        # MM/D/YYYY
        # DD/M/YYYY
        if re.search("^\d{2}(/|-)\d{1}(/|-)\d{4}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = raw_date[6:]
            month, day = check(month, day, raw_date)

        # ?/?/YYYY
        # M/D/YYYY
        # D/M/YYYY
        if re.search("^\d{1}(/|-)\d{1}(/|-)\d{4}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = raw_date[6:]
            month, day = check(month, day, raw_date)

        # ?/??/YY
        # M/DD/YY
        # D/MM/YY
        if re.search("^\d{1}(/|-)\d{2}(/|-)\d{2}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = raw_date[6:]
            month, day = check(month, day, raw_date)
            raise Exception(f"Unexpected date format: {raw_date}")

        # ??/?/YY
        # MM/D/YY
        # DD/M/YY
        if re.search("^\d{2}(/|-)\d{1}(/|-)\d{2}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = "19" + raw_date[6:]
            month, day = check(month, day, raw_date)
            raise Exception(f"Unexpected date format: {raw_date}")

        # ?/?/YY
        # M/D/YY
        # D/M/YY
        if re.search("^\d{1}(/|-)\d{1}(/|-)\d{2}$", raw_date):
            day = int(raw_date[:2])
            month = int(raw_date[3:5])
            year = "19" + raw_date[6:]
            month, day = check(month, day, raw_date)
            raise Exception(f"Unexpected date format: {raw_date}")


        # Two dates ################################################

        # YYYY/YYYY
        if re.search("^\d{4}(/|-)\d{4}$", raw_date):
            year = int(raw_date[:4])
            year_2 = int(raw_date[5:9])
            if year_2 < year:
                year = year_2

        # YYYY-MM/YYYY-MM
        search_results = re.search(
            "^(\d{4})(/|-)(\d{2})"
            "(/|-)"
            "(\d{4})(/|-)(\d{2})$",
            raw_date
        )
        if search_results:
            year = search_results.group(1)
            month = search_results.group(3)
            year_2 = search_results.group(5)
            month_2 = search_results.group(7)
            if int(year_2) + int(month_2) < int(year) + int(month):
                month = month_2
                year = year_2

        # YYYY-MM-DD/YYYY-MM-DD
        search_results = re.search(
            "^(\d{4})(/|-)(\d{2})(/|-)(\d{2})"
            "(/|-)"
            "(\d{4})(/|-)(\d{2})(/|-)(\d{2})$",
            raw_date
        )
        if search_results:
            year = search_results.group(1)
            month = search_results.group(3)
            day = search_results.group(5)
            year_2 = search_results.group(7)
            month_2 = search_results.group(9)
            day_2 = search_results.group(11)
            if int(year_2) + int(month_2) + int(day_2) < int(year) + int(month) + int(day):
                month = month_2
                year = year_2
                day = day_2

        # M...-YYYY/M...-YYYY
        search_results = re.search(
            "^(" + MONTHS_REGEX + ")(/|-)(\d{4})"
            "(/|-)(" + MONTHS_REGEX + ")(/|-)(\d{4})$"
            , raw_date
        )
        if search_results:
            month = MONTH_NAME_MAP[search_results.group(1).lower()]
            year = search_results.group(3)
            month_2 = MONTH_NAME_MAP[search_results.group(5).lower()]
            year_2 = search_results.group(7)
            if int(year_2) + int(month_2) < int(year) + int(month):
                month = month_2
                year = year_2

        # M...-YYYY/DD-M...-YYYY
        search_results = re.search(
            "^(" + MONTHS_REGEX + ")(/|-)(\d{4})"
            "(/|-)"
            "(\d{2})(/|-)(" + MONTHS_REGEX + ")(/|-)(\d{4})$"
            , raw_date
        )
        if search_results:
            month = MONTH_NAME_MAP[search_results.group(1).lower()]
            year = search_results.group(3)
            day = search_results.group(5)
            month_2 = MONTH_NAME_MAP[search_results.group(7).lower()]
            year_2 = search_results.group(9)
            if int(year_2) + int(month_2) < int(year) + int(month):
                month = month_2
                year = year_2

        # M... D, YYYY
        search_results = re.search(
            "^(" + MONTHS_REGEX + ") "
            "(\d|\d{2}), "
            "(\d{4})$"
            , raw_date
        )
        if search_results:
            month = MONTH_NAME_MAP[search_results.group(1).lower()]
            day = search_results.group(2)
            year = search_results.group(3)

        normalized_date = f"{year}-{month}-{day}"

        if normalized_date == "0000-00-00":
            raise Exception(f"Unexpected date format: {raw_date}\n")

        return normalized_date

