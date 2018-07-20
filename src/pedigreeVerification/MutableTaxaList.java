package pedigreeVerification;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;

public class MutableTaxaList implements TaxaList {
    private ArrayList<Taxon> myTaxaList;
    
    MutableTaxaList(){
        myTaxaList = new ArrayList<Taxon>();
    }
    
    @Override
    public int size() {
        // TODO Auto-generated method stub
        return myTaxaList.size();
    }

    @Override
    public boolean isEmpty() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean contains(Object o) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public Iterator<Taxon> iterator() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Object[] toArray() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public <T> T[] toArray(T[] a) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean add(Taxon e) {
        // TODO Auto-generated method stub
        myTaxaList.add(e);
        return true;
    }

    @Override
    public boolean remove(Object o) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean addAll(Collection<? extends Taxon> c) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean addAll(int index, Collection<? extends Taxon> c) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public void clear() {
        // TODO Auto-generated method stub

    }

    @Override
    public Taxon get(int index) {
        
        return myTaxaList.get(index);
    }

    @Override
    public Taxon set(int index, Taxon element) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void add(int index, Taxon element) {
        // TODO Auto-generated method stub

    }

    @Override
    public Taxon remove(int index) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int indexOf(Object o) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int lastIndexOf(Object o) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public ListIterator<Taxon> listIterator() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public ListIterator<Taxon> listIterator(int index) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<Taxon> subList(int fromIndex, int toIndex) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int numberOfTaxa() {
        // TODO Auto-generated method stub
        
        return this.myTaxaList.size();
    }

    @Override
    public String taxaName(int index) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int indexOf(String name) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int indexOf(Taxon taxon) {
        // TODO Auto-generated method stub
        return 0;
    }

}
