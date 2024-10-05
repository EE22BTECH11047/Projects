import CONFIG from '../config';

class Storage {
  static getItems(type = null, overwrite_items = null) {
    let items =
      overwrite_items !== null
        ? JSON.stringify(overwrite_items)
        : localStorage.getItem(CONFIG.LOCALSTORAGE_KEY);
    let filtered_items = [];

    try {
      items = JSON.parse(items);
    } catch (e) {
      items = [];
    }

    items = items || [];

    items.forEach(function (item) {
      item.amount = parseInt(item.amount);
      item.apr = parseInt(item.apr);

      if (type && item.type !== type) {
        return item;
      }

      filtered_items.push(item);

      return item;
    });

    return filtered_items;
  }

  static setItems(items) {
    localStorage.setItem(CONFIG.LOCALSTORAGE_KEY, JSON.stringify(items));
  }

  static addItem(item) {
    let items = Storage.getItems();

    items.push(item);

    Storage.setItems(items);
  }

  static removeItem(id) {
    let items = Storage.getItems();

    const index = items.findIndex((item) => item.id === id);

    if (index >= 0) {
      items.splice(index, 1);
    }

    Storage.setItems(items);
  }
}

export default Storage;
