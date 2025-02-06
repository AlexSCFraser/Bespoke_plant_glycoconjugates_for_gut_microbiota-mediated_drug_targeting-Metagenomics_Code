import re


def adjusted_labels(labels, name_title_insert=None):
    # todo: remove this renaming here, use adjusted_names instead and specify the angle in caller.
    for item in labels:
        if name_title_insert:
            item.set_text('\n'.join(item._text.split(' ')).replace(f"{name_title_insert}\n", ''))
        else:
            item.set_text('\n'.join(item._text.split(' ')))
        item.set_text(item._text.replace('\ncoverage', '').replace('\npresent', ''))

    if len(labels) > 9:
        pul_name_angle = 45
        # x_labels = labels
    else:
        pul_name_angle = 0

    return pul_name_angle, labels


def adjusted_names(names: list[str], name_title_insert=None):
    new_names = []
    for item in names:
            new_names.append(re.sub(f'(?i)([\n_]{name_title_insert})?[\n_](coverage|present)', '', item).capitalize())

    return new_names
